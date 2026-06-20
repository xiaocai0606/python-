#!/usr/bin/env python3
"""RINEX 2.x GPS pseudorange single point positioning.

The program reads a RINEX observation file (*.o) and navigation file (*.n),
then solves receiver ECEF coordinates and receiver clock bias by iterative
least squares using broadcast ephemerides.

本程序完成一个基础 GPS 伪距单点定位流程：
1. 读取 RINEX O 文件中的伪距观测值；
2. 读取 RINEX N 文件中的广播星历；
3. 根据广播星历计算卫星发射时刻的位置和钟差；
4. 对每个历元用最小二乘迭代求接收机 XYZ 坐标和接收机钟差。
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# GPS 常用常数。程序里所有距离单位为 m，时间单位为 s，角度计算时用 rad。
C = 299_792_458.0  # 真空光速
MU = 3.986005e14  # WGS-84 地球引力常数
OMEGA_E = 7.2921151467e-5  # 地球自转角速度
F_REL = -4.442807633e-10  # 广播星历相对论钟差改正常数
GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)  # GPS 周秒起算历元


@dataclass
class NavHeader:
    """N 文件头部信息，只保留定位中可能用到的几项。"""

    version: str = ""
    ion_alpha: Tuple[float, ...] = ()
    ion_beta: Tuple[float, ...] = ()
    leap_seconds: int = 0


@dataclass
class Ephemeris:
    """一颗 GPS 卫星的一组广播星历参数。

    RINEX 2.x N 文件中每颗卫星一组星历占 8 行，这里的字段基本按
    广播星历顺序保存。后续会用这些参数计算卫星位置和卫星钟差。
    """

    prn: int
    epoch: datetime
    toc: float
    af0: float
    af1: float
    af2: float
    iode: float
    crs: float
    delta_n: float
    m0: float
    cuc: float
    ecc: float
    cus: float
    sqrt_a: float
    toe: float
    cic: float
    omega0: float
    cis: float
    i0: float
    crc: float
    arg_perigee: float
    omega_dot: float
    idot: float
    gps_week: int
    tgd: float


@dataclass
class ObsHeader:
    """O 文件头部信息。

    approx_position 是文件头给出的测站概略坐标，可作为迭代初值。
    obs_types 是观测类型列表，例如 L1、L2、C1、P1、P2、S1、S2。
    """

    version: str = ""
    marker_name: str = ""
    approx_position: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    obs_types: List[str] = None
    interval: float = 0.0
    first_obs: Optional[datetime] = None

    def __post_init__(self) -> None:
        if self.obs_types is None:
            self.obs_types = []


@dataclass
class ObservationEpoch:
    """一个观测历元的数据。

    observations 的结构为：
        {"G21": {"C1": 22064095.016, "P2": 22064104.152, ...}, ...}
    """

    epoch: datetime
    gps_week: int
    sow: float
    flag: int
    observations: Dict[str, Dict[str, float]]


@dataclass
class SatelliteUse:
    """记录某次定位中实际参与解算的一颗卫星及其残差。"""

    sat: str
    pseudorange: float
    x: float
    y: float
    z: float
    clock_bias: float
    residual: float
    elevation_deg: Optional[float] = None


@dataclass
class SppSolution:
    """一个历元的最终定位结果。"""

    epoch: datetime
    x: float
    y: float
    z: float
    receiver_clock_s: float
    iterations: int
    used: List[SatelliteUse]
    rms: float
    lat_deg: float
    lon_deg: float
    height_m: float


def parse_float(text: str, default: float = 0.0) -> float:
    """解析 RINEX 中的浮点数。

    RINEX 常用 D 表示科学计数法指数，例如 0.123D-03，Python 需要
    先把 D 换成 E 才能转为 float。
    """

    text = text.strip().replace("D", "E").replace("d", "E")
    if not text:
        return default
    return float(text)


def parse_year(two_or_four_digit: int) -> int:
    """把 RINEX 2.x 中的两位年份转成四位年份。"""

    if two_or_four_digit < 80:
        return 2000 + two_or_four_digit
    if two_or_four_digit < 100:
        return 1900 + two_or_four_digit
    return two_or_four_digit


def datetime_to_gps(dt: datetime) -> Tuple[int, float]:
    """把 UTC datetime 转成 GPS 周和周内秒。

    本实验数据为 GPS 时间系统，未在这里额外处理闰秒。
    """

    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    delta = dt - GPS_EPOCH
    total = delta.total_seconds()
    week = int(total // 604800)
    return week, total - week * 604800


def timediff_seconds(t: float) -> float:
    """处理 GPS 周跳变，使时间差落在半周范围内。"""

    if t > 302400.0:
        t -= 604800.0
    elif t < -302400.0:
        t += 604800.0
    return t


def read_nav(path: Path) -> Tuple[NavHeader, Dict[int, List[Ephemeris]]]:
    """读取 RINEX 2.x 导航文件。

    返回值中第二项按 PRN 编号分组，例如 nav[21] 表示 G21 的全部星历。
    同一颗卫星一天内会有多组星历，定位时再按 TOE 最近原则选择。
    """

    header = NavHeader()
    eph_by_prn: Dict[int, List[Ephemeris]] = {}

    with path.open("r", encoding="ascii", errors="ignore") as f:
        # 先读取文件头。RINEX 文件头每行标签通常在第 61 列之后。
        for line in f:
            label = line[60:].strip() if len(line) >= 60 else ""
            if "RINEX VERSION / TYPE" in label:
                header.version = line[:20].strip()
            elif "ION ALPHA" in label:
                header.ion_alpha = tuple(parse_float(item) for item in line[:60].split()[:4])
            elif "ION BETA" in label:
                header.ion_beta = tuple(parse_float(item) for item in line[:60].split()[:4])
            elif "LEAP SECONDS" in label:
                header.leap_seconds = int(parse_float(line[:6]))
            elif "END OF HEADER" in label:
                break

        # 文件头后面就是星历记录。GPS RINEX 2.x 中每组星历固定 8 行：
        # 第 1 行为 PRN、钟参考时刻和卫星钟参数，后 7 行为轨道参数。
        while True:
            first = f.readline()
            if not first:
                break
            if not first.strip():
                continue
            block = [first] + [f.readline() for _ in range(7)]
            if len(block[-1]) == 0:
                break
            try:
                # 第 1 行：卫星号、星历参考历元、卫星钟差多项式参数。
                prn = int(block[0][0:2])
                year = parse_year(int(block[0][2:5]))
                month = int(block[0][5:8])
                day = int(block[0][8:11])
                hour = int(block[0][11:14])
                minute = int(block[0][14:17])
                second = parse_float(block[0][17:22])
                sec_int = int(second)
                micro = int(round((second - sec_int) * 1_000_000))
                epoch = datetime(year, month, day, hour, minute, sec_int, micro, tzinfo=timezone.utc)
                _, toc = datetime_to_gps(epoch)

                # 星历主体按固定列宽 19 字符读取。第一行从第 23 列开始，
                # 后 7 行从第 4 列开始，每行最多 4 个参数。
                values = [
                    [parse_float(block[0][22:41]), parse_float(block[0][41:60]), parse_float(block[0][60:79])],
                ]
                for line in block[1:]:
                    values.append([parse_float(line[i : i + 19]) for i in range(3, 79, 19)])

                eph = Ephemeris(
                    prn=prn,
                    epoch=epoch,
                    toc=toc,
                    af0=values[0][0],
                    af1=values[0][1],
                    af2=values[0][2],
                    iode=values[1][0],
                    crs=values[1][1],
                    delta_n=values[1][2],
                    m0=values[1][3],
                    cuc=values[2][0],
                    ecc=values[2][1],
                    cus=values[2][2],
                    sqrt_a=values[2][3],
                    toe=values[3][0],
                    cic=values[3][1],
                    omega0=values[3][2],
                    cis=values[3][3],
                    i0=values[4][0],
                    crc=values[4][1],
                    arg_perigee=values[4][2],
                    omega_dot=values[4][3],
                    idot=values[5][0],
                    gps_week=int(values[5][2]),
                    tgd=values[6][2],
                )
            except (ValueError, IndexError):
                # 若遇到不完整或格式异常记录，跳过该组星历。
                continue
            eph_by_prn.setdefault(prn, []).append(eph)

    # 每颗卫星内部按 TOE 排序，便于后续选取最近星历。
    for records in eph_by_prn.values():
        records.sort(key=lambda e: e.toe)
    return header, eph_by_prn


def read_obs(path: Path) -> Tuple[ObsHeader, List[ObservationEpoch]]:
    """读取 RINEX 2.x 观测文件。

    本程序主要使用伪距观测值，例如 C1、P1、P2。载波相位 L1/L2
    和信噪比 S1/S2 会读入，但不参与基础伪距单点定位。
    """

    header = ObsHeader()
    epochs: List[ObservationEpoch] = []

    with path.open("r", encoding="ascii", errors="ignore") as f:
        # 读取 O 文件头：测站名、概略坐标、观测类型、采样间隔等。
        for line in f:
            label = line[60:].strip() if len(line) >= 60 else ""
            if "RINEX VERSION / TYPE" in label:
                header.version = line[:20].strip()
            elif "MARKER NAME" in label:
                header.marker_name = line[:60].strip()
            elif "APPROX POSITION XYZ" in label:
                header.approx_position = (
                    parse_float(line[0:14]),
                    parse_float(line[14:28]),
                    parse_float(line[28:42]),
                )
            elif "# / TYPES OF OBSERV" in label:
                count = int(parse_float(line[:6]))
                obs_types = line[6:60].split()
                # 观测类型较多时，RINEX 会把类型列表续写到下一行。
                while len(obs_types) < count:
                    cont = f.readline()
                    obs_types.extend(cont[6:60].split())
                header.obs_types = obs_types[:count]
            elif "INTERVAL" in label:
                header.interval = parse_float(line[:10])
            elif "TIME OF FIRST OBS" in label:
                year = int(parse_float(line[0:6]))
                month = int(parse_float(line[6:12]))
                day = int(parse_float(line[12:18]))
                hour = int(parse_float(line[18:24]))
                minute = int(parse_float(line[24:30]))
                second = parse_float(line[30:43])
                sec_int = int(second)
                micro = int(round((second - sec_int) * 1_000_000))
                header.first_obs = datetime(year, month, day, hour, minute, sec_int, micro, tzinfo=timezone.utc)
            elif "END OF HEADER" in label:
                break

        # 文件头之后按历元读取。每个历元先有一行历元头：
        # 年月日时分秒、历元标志、卫星数、卫星编号列表。
        while True:
            line = f.readline()
            if not line:
                break
            if len(line) < 32 or not line[:3].strip():
                continue
            try:
                year = parse_year(int(line[1:3]))
                month = int(line[4:6])
                day = int(line[7:9])
                hour = int(line[10:12])
                minute = int(line[13:15])
                second = parse_float(line[16:26])
                flag = int(line[28:29])
                nsat = int(line[29:32])
            except ValueError:
                continue

            # 如果卫星数量很多，卫星编号列表可能跨多行。
            sat_text = line[32:].rstrip("\n")
            while len(sat_text.replace(" ", "")) < nsat * 3:
                sat_text += f.readline().rstrip("\n")
            sats = [sat_text[i : i + 3].strip() for i in range(0, nsat * 3, 3)]

            sec_int = int(second)
            micro = int(round((second - sec_int) * 1_000_000))
            epoch_dt = datetime(year, month, day, hour, minute, sec_int, micro, tzinfo=timezone.utc)
            gps_week, sow = datetime_to_gps(epoch_dt)

            # flag 为 0 或 1 时一般表示正常观测；其他类型可能是事件记录等，
            # 这里跳过并读掉对应观测行，避免后续错位。
            if flag not in (0, 1):
                for _ in sats:
                    for _ in range((len(header.obs_types) + 4) // 5):
                        f.readline()
                continue

            observations: Dict[str, Dict[str, float]] = {}
            lines_per_sat = (len(header.obs_types) + 4) // 5
            for sat in sats:
                chunks = ""
                # 每颗卫星每行最多 5 个观测值，每个观测值占 16 字符。
                for _ in range(lines_per_sat):
                    chunks += f.readline().rstrip("\n").ljust(80)
                sat_obs: Dict[str, float] = {}
                for idx, obs_type in enumerate(header.obs_types):
                    field = chunks[idx * 16 : (idx + 1) * 16]
                    value = field[:14].strip()
                    if value:
                        sat_obs[obs_type] = float(value)
                observations[sat] = sat_obs

            epochs.append(ObservationEpoch(epoch_dt, gps_week, sow, flag, observations))

    return header, epochs


def choose_ephemeris(records: Dict[int, List[Ephemeris]], prn: int, transmit_sow: float) -> Optional[Ephemeris]:
    """为某颗卫星选择 TOE 离信号发射时刻最近的一组星历。"""

    candidates = records.get(prn, [])
    if not candidates:
        return None
    return min(candidates, key=lambda e: abs(timediff_seconds(transmit_sow - e.toe)))


def kepler_eccentric_anomaly(mean_anomaly: float, ecc: float) -> float:
    """用迭代法解 Kepler 方程，得到偏近点角 E。"""

    eccentric = mean_anomaly
    for _ in range(30):
        nxt = mean_anomaly + ecc * math.sin(eccentric)
        if abs(nxt - eccentric) < 1e-13:
            return nxt
        eccentric = nxt
    return eccentric


def satellite_clock(eph: Ephemeris, transmit_sow: float, eccentric_anomaly: Optional[float] = None) -> float:
    """计算卫星钟差。

    包含广播星历钟差多项式、相对论效应改正和 TGD 改正。
    返回值单位为秒。
    """

    dt = timediff_seconds(transmit_sow - eph.toc)
    clock = eph.af0 + eph.af1 * dt + eph.af2 * dt * dt
    if eccentric_anomaly is not None:
        clock += F_REL * eph.ecc * eph.sqrt_a * math.sin(eccentric_anomaly)
    clock -= eph.tgd
    return clock


def satellite_position(eph: Ephemeris, transmit_sow: float) -> Tuple[float, float, float, float]:
    """根据广播星历计算卫星在发射时刻的 ECEF 坐标和钟差。

    这部分对应 GPS 广播星历标准计算流程：平均角速度 -> 平近点角 ->
    偏近点角 -> 真近点角 -> 摄动改正 -> 轨道平面坐标 -> 地固坐标。
    """

    a = eph.sqrt_a * eph.sqrt_a
    n0 = math.sqrt(MU / (a * a * a))
    tk = timediff_seconds(transmit_sow - eph.toe)
    n = n0 + eph.delta_n
    mean_anomaly = eph.m0 + n * tk
    eccentric = kepler_eccentric_anomaly(mean_anomaly, eph.ecc)

    # 由偏近点角 E 计算真近点角 v。
    sin_v = math.sqrt(1.0 - eph.ecc * eph.ecc) * math.sin(eccentric) / (1.0 - eph.ecc * math.cos(eccentric))
    cos_v = (math.cos(eccentric) - eph.ecc) / (1.0 - eph.ecc * math.cos(eccentric))
    true_anomaly = math.atan2(sin_v, cos_v)
    phi = true_anomaly + eph.arg_perigee

    # 二阶调和项改正：改正升交距角、轨道半径和轨道倾角。
    du = eph.cus * math.sin(2.0 * phi) + eph.cuc * math.cos(2.0 * phi)
    dr = eph.crs * math.sin(2.0 * phi) + eph.crc * math.cos(2.0 * phi)
    di = eph.cis * math.sin(2.0 * phi) + eph.cic * math.cos(2.0 * phi)

    u = phi + du
    r = a * (1.0 - eph.ecc * math.cos(eccentric)) + dr
    inc = eph.i0 + di + eph.idot * tk

    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)
    # 升交点赤经需要扣除地球自转影响。
    omega = eph.omega0 + (eph.omega_dot - OMEGA_E) * tk - OMEGA_E * eph.toe

    # 从轨道平面坐标转换到地心地固坐标 ECEF。
    x = x_orb * math.cos(omega) - y_orb * math.cos(inc) * math.sin(omega)
    y = x_orb * math.sin(omega) + y_orb * math.cos(inc) * math.cos(omega)
    z = y_orb * math.sin(inc)
    return x, y, z, satellite_clock(eph, transmit_sow, eccentric)


def rotate_for_earth(x: float, y: float, z: float, travel_time: float) -> Tuple[float, float, float]:
    """地球自转改正。

    信号从卫星传播到接收机需要约 0.07 s，这段时间地球仍在旋转。
    因此要把发射时刻的卫星坐标旋转到接收时刻对应的地固坐标系。
    """

    angle = OMEGA_E * travel_time
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    return cos_a * x + sin_a * y, -sin_a * x + cos_a * y, z


def solve_linear_4(normal: List[List[float]], rhs: List[float]) -> List[float]:
    """解 4 阶线性方程组。

    最小二乘法最后会得到 4 个未知数的法方程：
        normal * dx = rhs
    这里用高斯消元求解 dx = [dX, dY, dZ, dClock_m]。
    """

    a = [row[:] + [rhs[i]] for i, row in enumerate(normal)]
    n = 4
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1e-18:
            raise ValueError("normal equation is singular")
        a[col], a[pivot] = a[pivot], a[col]
        div = a[col][col]
        for j in range(col, n + 1):
            a[col][j] /= div
        for r in range(n):
            if r == col:
                continue
            factor = a[r][col]
            for j in range(col, n + 1):
                a[r][j] -= factor * a[col][j]
    return [a[i][n] for i in range(n)]


def least_squares(rows: List[List[float]], obs_minus_calc: List[float]) -> List[float]:
    """由误差方程组成法方程并求最小二乘解。

    rows 是设计矩阵 B 的每一行，obs_minus_calc 是 O-C 向量 L。
    解算公式等价于：
        dx = (B^T B)^-1 B^T L
    """

    normal = [[0.0] * 4 for _ in range(4)]
    rhs = [0.0] * 4
    for row, value in zip(rows, obs_minus_calc):
        for i in range(4):
            rhs[i] += row[i] * value
            for j in range(4):
                normal[i][j] += row[i] * row[j]
    return solve_linear_4(normal, rhs)


def ecef_to_geodetic(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """WGS-84 地心地固坐标 XYZ 转大地坐标 BLH。

    返回纬度、经度和大地高。纬度和经度单位为度，高程单位为 m。
    """

    a = 6378137.0
    f = 1.0 / 298.257223563
    e2 = f * (2.0 - f)
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1.0 - e2))
    h = 0.0
    for _ in range(20):
        sin_lat = math.sin(lat)
        n = a / math.sqrt(1.0 - e2 * sin_lat * sin_lat)
        h = p / math.cos(lat) - n
        nxt = math.atan2(z, p * (1.0 - e2 * n / (n + h)))
        if abs(nxt - lat) < 1e-13:
            lat = nxt
            break
        lat = nxt
    return math.degrees(lat), math.degrees(lon), h


def elevation_angle(receiver: Sequence[float], satellite: Sequence[float]) -> float:
    """计算卫星相对接收机的高度角，主要用于可选的高度角截止。"""

    lat, lon, _ = ecef_to_geodetic(receiver[0], receiver[1], receiver[2])
    lat_r = math.radians(lat)
    lon_r = math.radians(lon)
    dx = satellite[0] - receiver[0]
    dy = satellite[1] - receiver[1]
    dz = satellite[2] - receiver[2]
    east = -math.sin(lon_r) * dx + math.cos(lon_r) * dy
    north = -math.sin(lat_r) * math.cos(lon_r) * dx - math.sin(lat_r) * math.sin(lon_r) * dy + math.cos(lat_r) * dz
    up = math.cos(lat_r) * math.cos(lon_r) * dx + math.cos(lat_r) * math.sin(lon_r) * dy + math.sin(lat_r) * dz
    return math.degrees(math.atan2(up, math.hypot(east, north)))


def select_pseudorange(obs: Dict[str, float], preferred: Sequence[str]) -> Optional[float]:
    """按优先级选择一个可用伪距观测值。

    默认优先级是 C1、P1、P2。实际文件中某些观测类型可能为空，
    所以需要按顺序寻找第一个有效值。
    """

    for key in preferred:
        value = obs.get(key)
        if value and value > 1.0:
            return value
    return None


def solve_epoch(
    epoch: ObservationEpoch,
    nav: Dict[int, List[Ephemeris]],
    initial_xyz: Tuple[float, float, float],
    preferred_obs: Sequence[str],
    min_elevation: Optional[float] = None,
    max_iterations: int = 10,
) -> Optional[SppSolution]:
    """解算单个历元的伪距单点定位结果。

    未知数为接收机坐标 X/Y/Z 和接收机钟差。程序内部把接收机钟差
    先用距离单位 clock_m = c * dt 表示，这样可以和伪距残差同单位。
    """

    x, y, z = initial_xyz
    clock_m = 0.0
    used_final: List[SatelliteUse] = []

    # 外层循环是非线性最小二乘迭代。每次用当前近似坐标重新计算
    # 卫星几何距离和设计矩阵，然后求坐标改正数。
    for iteration in range(1, max_iterations + 1):
        rows: List[List[float]] = []
        residuals: List[float] = []
        used: List[SatelliteUse] = []

        for sat, obs_values in epoch.observations.items():
            # 本程序只处理 GPS 卫星，RINEX 中 GPS 卫星编号形如 G21。
            if not sat.startswith("G"):
                continue
            try:
                prn = int(sat[1:])
            except ValueError:
                continue
            pseudorange = select_pseudorange(obs_values, preferred_obs)
            if pseudorange is None:
                continue

            # 用伪距和当前接收机钟差估计信号发射时刻。
            # 第一次估计用于选星历和粗算卫星钟差。
            transmit_sow = epoch.sow - (pseudorange - clock_m) / C
            eph = choose_ephemeris(nav, prn, transmit_sow)
            if eph is None:
                continue
            sat_clock = satellite_clock(eph, transmit_sow)
            # 用卫星钟差修正发射时刻，然后再精算卫星坐标。
            transmit_sow = epoch.sow - (pseudorange - clock_m) / C - sat_clock

            sx, sy, sz, sat_clock = satellite_position(eph, transmit_sow)
            travel = (pseudorange - clock_m + C * sat_clock) / C
            sx, sy, sz = rotate_for_earth(sx, sy, sz, travel)
            rho = math.sqrt((sx - x) ** 2 + (sy - y) ** 2 + (sz - z) ** 2)
            if rho <= 0.0:
                continue

            elev = None
            if min_elevation is not None and any(abs(v) > 1.0 for v in (x, y, z)):
                elev = elevation_angle((x, y, z), (sx, sy, sz))
                if elev < min_elevation:
                    continue

            # 设计矩阵一行：几何距离对 X/Y/Z 的偏导数，以及钟差项。
            # 这里钟差未知数用 m 表示，所以最后一列系数为 1。
            row = [(x - sx) / rho, (y - sy) / rho, (z - sz) / rho, 1.0]
            # O-C：伪距观测值 - 计算距离 + 卫星钟差改正 - 当前接收机钟差。
            omc = pseudorange - rho + C * sat_clock - clock_m
            rows.append(row)
            residuals.append(omc)
            used.append(SatelliteUse(sat, pseudorange, sx, sy, sz, sat_clock, omc, elev))

        if len(rows) < 4:
            # 四个未知数至少需要四颗卫星。
            return None

        # 求本次迭代的坐标改正和钟差改正。
        dx, dy, dz, d_clock_m = least_squares(rows, residuals)
        x += dx
        y += dy
        z += dz
        clock_m += d_clock_m
        used_final = used

        # 坐标和钟差改正都足够小时认为迭代收敛。
        if math.sqrt(dx * dx + dy * dy + dz * dz) < 1e-4 and abs(d_clock_m) < 1e-4:
            break

    # 用最终坐标重新计算残差和 RMS，便于评价该历元定位质量。
    postfit: List[float] = []
    refreshed: List[SatelliteUse] = []
    for item in used_final:
        rho = math.sqrt((item.x - x) ** 2 + (item.y - y) ** 2 + (item.z - z) ** 2)
        res = item.pseudorange - rho + C * item.clock_bias - clock_m
        elev = elevation_angle((x, y, z), (item.x, item.y, item.z))
        refreshed.append(SatelliteUse(item.sat, item.pseudorange, item.x, item.y, item.z, item.clock_bias, res, elev))
        postfit.append(res)

    dof = max(1, len(postfit) - 4)
    rms = math.sqrt(sum(r * r for r in postfit) / dof)
    lat, lon, height = ecef_to_geodetic(x, y, z)
    return SppSolution(epoch.epoch, x, y, z, clock_m / C, iteration, refreshed, rms, lat, lon, height)


def write_csv(path: Path, solutions: Iterable[SppSolution]) -> None:
    """把一个或多个历元的定位结果保存为 CSV。"""

    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "epoch",
            "x_m",
            "y_m",
            "z_m",
            "lat_deg",
            "lon_deg",
            "height_m",
            "receiver_clock_s",
            "satellites",
            "rms_m",
            "iterations",
        ])
        for sol in solutions:
            writer.writerow([
                sol.epoch.isoformat(),
                f"{sol.x:.4f}",
                f"{sol.y:.4f}",
                f"{sol.z:.4f}",
                f"{sol.lat_deg:.10f}",
                f"{sol.lon_deg:.10f}",
                f"{sol.height_m:.4f}",
                f"{sol.receiver_clock_s:.12e}",
                len(sol.used),
                f"{sol.rms:.4f}",
                sol.iterations,
            ])


def main() -> int:
    """命令行入口：读取参数、读取文件、执行定位并输出结果。"""

    parser = argparse.ArgumentParser(description="RINEX GPS pseudorange single point positioning")
    parser.add_argument("obs_file", type=Path, help="RINEX observation file, for example wuhn1660.12o")
    parser.add_argument("nav_file", type=Path, help="RINEX navigation file, for example brdc1660.12n")
    parser.add_argument("--epoch", type=int, default=0, help="zero-based epoch index to solve")
    parser.add_argument("--all", action="store_true", help="solve all available epochs")
    parser.add_argument("--max-epochs", type=int, default=0, help="limit epochs when --all is used")
    parser.add_argument("--obs", default="C1,P1,P2", help="comma-separated pseudorange priority")
    parser.add_argument("--initial", choices=("approx", "zero"), default="approx", help="initial receiver coordinates")
    parser.add_argument("--min-elevation", type=float, default=None, help="optional elevation mask in degrees")
    parser.add_argument("--csv", type=Path, default=None, help="write solutions to CSV")
    args = parser.parse_args()

    obs_header, epochs = read_obs(args.obs_file)
    _, nav = read_nav(args.nav_file)
    preferred = [item.strip() for item in args.obs.split(",") if item.strip()]
    # 默认用 O 文件头中的概略坐标作为初值；也可以通过 --initial zero
    # 从 (0, 0, 0) 开始迭代，和教学流程中的初值设为 0 对应。
    initial = obs_header.approx_position if args.initial == "approx" else (0.0, 0.0, 0.0)

    print(f"Observation file: {args.obs_file}")
    print(f"Navigation file:  {args.nav_file}")
    print(f"Marker: {obs_header.marker_name or '(unknown)'}")
    print(f"Observation types: {' '.join(obs_header.obs_types)}")
    print(f"Epochs read: {len(epochs)}, navigation satellites: {len(nav)}")
    print(f"Initial XYZ: {initial[0]:.4f}, {initial[1]:.4f}, {initial[2]:.4f}")

    target_epochs = epochs if args.all else epochs[args.epoch : args.epoch + 1]
    if args.all and args.max_epochs > 0:
        target_epochs = target_epochs[: args.max_epochs]

    # 逐历元独立定位。单点定位不依赖前后历元，一个历元卫星数足够即可解。
    solutions: List[SppSolution] = []
    for epoch in target_epochs:
        solution = solve_epoch(epoch, nav, initial, preferred, args.min_elevation)
        if solution is not None:
            solutions.append(solution)

    if not solutions:
        print("No solution: fewer than four usable GPS pseudorange observations or missing ephemerides.")
        return 2

    for sol in solutions[:10]:
        print()
        print(f"Epoch: {sol.epoch.isoformat()}")
        print(f"XYZ (m): {sol.x:.4f}, {sol.y:.4f}, {sol.z:.4f}")
        print(f"BLH: lat {sol.lat_deg:.10f} deg, lon {sol.lon_deg:.10f} deg, h {sol.height_m:.4f} m")
        print(f"Receiver clock: {sol.receiver_clock_s:.12e} s")
        print(f"Used satellites: {len(sol.used)}, iterations: {sol.iterations}, RMS: {sol.rms:.4f} m")
        print("Satellites:", " ".join(item.sat for item in sol.used))

    if len(solutions) > 10:
        print(f"\n... {len(solutions) - 10} more solutions omitted from screen output")

    if args.csv:
        write_csv(args.csv, solutions)
        print(f"\nCSV written: {args.csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
