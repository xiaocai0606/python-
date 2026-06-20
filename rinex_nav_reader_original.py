class RinexHeader:  # 用于存储 RINEX 文件的文件头信息
    def __init__(self):
        self.version = ""  # RINEX版本
        self.file_type = ""  # 文件类型
        self.pgm = ""  # 程序名称
        self.run_by = ""  # 运行机构
        self.date = ""  # 运行日期
        self.ion_alpha = []  # 电离层参数 alpha
        self.ion_beta = []  # 电离层参数 beta
        self.delta_utc = []  # UTC参数
        self.leap_seconds = 0  # 跳秒数


class NavData:  # 用于存储每颗卫星的导航数据
    def __init__(self):
        self.prn = 0  # 卫星PRN号
        self.epoch = []  # 历元时间
        self.sv_clock = []  # 卫星钟差参数
        self.broadcast_orbit_1 = []  # 广播轨道1
        self.broadcast_orbit_2 = []  # 广播轨道2
        self.broadcast_orbit_3 = []  # 广播轨道3
        self.broadcast_orbit_4 = []  # 广播轨道4
        self.broadcast_orbit_5 = []  # 广播轨道5
        self.broadcast_orbit_6 = []  # 广播轨道6
        self.broadcast_orbit_7 = []  # 广播轨道7


class RinexNavReader:  # 读取和解析 RINEX 导航文件
    def __init__(self):
        self.header = RinexHeader()
        self.nav_data = []

    def _parse_scientific(self, val_str):
        """解析科学计数法格式的数值"""
        try:
            val_str = val_str.strip().replace('D', 'E')
            return float(val_str)
        except ValueError:
            return 0.0

    def read_file(self, filename):
        with open(filename, 'r') as f:
            # 读取文件头
            while True:
                line = f.readline()
                if "END OF HEADER" in line:
                    break
                self._parse_header_line(line)

            # 读取导航数据
            while True:
                nav_record = self._read_nav_record(f)
                if not nav_record:
                    break
                self.nav_data.append(nav_record)

    def _parse_header_line(self, line):
        if "RINEX VERSION / TYPE" in line:
            self.header.version = line[0:20].strip()
            self.header.file_type = line[20:40].strip()
        elif "PGM / RUN BY / DATE" in line:
            self.header.pgm = line[0:20].strip()
            self.header.run_by = line[20:40].strip()
            self.header.date = line[40:60].strip()
        elif "ION ALPHA" in line:
            self.header.ion_alpha = [self._parse_scientific(line[i:i + 12]) for i in range(0, 48, 12)]
        elif "ION BETA" in line:
            self.header.ion_beta = [self._parse_scientific(line[i:i + 12]) for i in range(0, 48, 12)]
        elif "DELTA-UTC" in line:
            self.header.delta_utc = [self._parse_scientific(line[i:i + 19]) for i in range(0, 57, 19)]
        elif "LEAP SECONDS" in line:
            self.header.leap_seconds = int(line[0:6])

    def _read_nav_record(self, f):
        # 读取第一行
        first_line = f.readline()
        if not first_line:
            return None

        nav = NavData()
        try:
            # 解析PRN号和历元
            nav.prn = int(first_line[0:2])
            nav.epoch = [
                int(first_line[2:5]),  # 年
                int(first_line[5:8]),  # 月
                int(first_line[8:11]),  # 日
                int(first_line[11:14]),  # 时
                int(first_line[14:17]),  # 分
                float(first_line[17:22])  # 秒
            ]

            # 解析卫星钟差参数
            nav.sv_clock = [
                self._parse_scientific(first_line[22:41]),  # 卫星钟偏差(s)
                self._parse_scientific(first_line[41:60]),  # 卫星钟漂移(s/s)
                self._parse_scientific(first_line[60:79])  # 卫星钟漂移速率(s/s^2)
            ]

            # 读取并解析后续7行轨道参数
            orbit_lines = []
            for _ in range(7):
                line = f.readline()
                if not line or len(line.strip()) < 4:
                    return None
                orbit_lines.append(line)

            # 解析轨道参数
            nav.broadcast_orbit_1 = [self._parse_scientific(orbit_lines[0][i:i + 19]) for i in range(3, 79, 19)]
            nav.broadcast_orbit_2 = [self._parse_scientific(orbit_lines[1][i:i + 19]) for i in range(3, 79, 19)]
            nav.broadcast_orbit_3 = [self._parse_scientific(orbit_lines[2][i:i + 19]) for i in range(3, 79, 19)]
            nav.broadcast_orbit_4 = [self._parse_scientific(orbit_lines[3][i:i + 19]) for i in range(3, 79, 19)]
            nav.broadcast_orbit_5 = [self._parse_scientific(orbit_lines[4][i:i + 19]) for i in range(3, 79, 19)]
            nav.broadcast_orbit_6 = [self._parse_scientific(orbit_lines[5][i:i + 19]) for i in range(3, 79, 19)]
            nav.broadcast_orbit_7 = [self._parse_scientific(orbit_lines[6][i:i + 19]) for i in range(3, 79, 19)]

        except Exception as e:
            print(f"警告：解析导航数据记录时出错: {e}")
            return None

        return nav

    def save_to_file(self, filename, csv_filename=None):
        """将数据保存到文件
        Args:
            filename: 保存文件头信息的文件名
            csv_filename: 保存导航数据的CSV文件名，如果为None则不保存CSV
        """
        # 保存文件头信息
        with open(filename, 'w', encoding='utf-8') as f:
            f.write("=== 文件头信息 ===\n")
            f.write(f"RINEX版本: {self.header.version}\n")
            f.write(f"文件类型: {self.header.file_type}\n")
            f.write(f"程序名称: {self.header.pgm}\n")
            f.write(f"运行机构: {self.header.run_by}\n")
            f.write(f"运行日期: {self.header.date}\n")
            f.write(f"电离层参数 alpha: {self.header.ion_alpha}\n")
            f.write(f"电离层参数 beta: {self.header.ion_beta}\n")
            f.write(f"UTC参数: {self.header.delta_utc}\n")
            f.write(f"跳秒数: {self.header.leap_seconds}\n")

        # 保存导航数据到CSV文件
        if csv_filename:
            import csv
            with open(csv_filename, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                # 写入CSV表头
                header = ['PRN', '历元', '钟偏差(s)', '钟漂移(s/s)', '钟漂移速率(s/s^2)',
                          '轨道参数1', '轨道参数2', '轨道参数3', '轨道参数4',
                          '轨道参数5', '轨道参数6', '轨道参数7']
                writer.writerow(header)

                # 写入导航数据
                for nav in self.nav_data:
                    epoch_str = f"{nav.epoch[0]:02d}/{nav.epoch[1]:02d}/{nav.epoch[2]:02d} " \
                                f"{nav.epoch[3]:02d}:{nav.epoch[4]:02d}:{nav.epoch[5]:06.3f}"

                    row = [
                        nav.prn,
                        epoch_str,
                        nav.sv_clock[0],
                        nav.sv_clock[1],
                        nav.sv_clock[2],
                        ','.join(map(str, nav.broadcast_orbit_1)),
                        ','.join(map(str, nav.broadcast_orbit_2)),
                        ','.join(map(str, nav.broadcast_orbit_3)),
                        ','.join(map(str, nav.broadcast_orbit_4)),
                        ','.join(map(str, nav.broadcast_orbit_5)),
                        ','.join(map(str, nav.broadcast_orbit_6)),
                        ','.join(map(str, nav.broadcast_orbit_7))
                    ]
                    writer.writerow(row)


# 使用示例
if __name__ == "__main__":
    reader = RinexNavReader()
    reader.read_file(r"C:\Users\xiaocai\Desktop\读取\brdc1660.12n")  # 替换为你的导航文件名

    # 保存数据到文件
    reader.save_to_file("导航文件头信息.txt", "导航数据.csv")