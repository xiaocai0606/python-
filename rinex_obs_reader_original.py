class RinexHeader:
    def __init__(self):
        self.version = ""  # RINEX版本
        self.file_type = ""  # 文件类型
        self.marker_name = ""  # 测站名
        self.marker_number = ""  # 测站编号
        self.observer = ""  # 观测者
        self.agency = ""  # 机构
        self.rec_number = ""  # 接收机编号
        self.rec_type = ""  # 接收机类型
        self.rec_version = ""  # 接收机版本
        self.ant_number = ""  # 天线编号
        self.ant_type = ""  # 天线类型
        self.approx_pos = []  # 概略坐标
        self.antenna_delta = []  # 天线高度
        self.obs_types = []  # 观测值类型
        self.interval = 0.0  # 观测间隔
        self.time_first_obs = []  # 首次观测时间
        self.leap_seconds = 0  # 跳秒数


class ObsData:
    def __init__(self):
        self.epoch = []  # 历元时间
        self.flag = 0  # 历元标志
        self.num_sats = 0  # 卫星数
        self.sats = []  # 卫星列表
        self.observations = {}  # 观测值字典 {PRN: {类型: {'value': 值, 'lli': LLI值, 'strength': 信号强度}}}


class RinexObsReader:
    def __init__(self):
        self.header = RinexHeader()
        self.obs_data = []

    def read_file(self, filename):
        with open(filename, 'r') as f:
            # 读取文件头
            while True:
                line = f.readline()
                if "END OF HEADER" in line:
                    break
                # 跳过COMMENT行
                if "COMMENT" in line.strip():
                    continue
                self._parse_header_line(line)

            # 读取观测数据
            while True:
                # 读取下一行
                line = f.readline()
                if not line:
                    break

                # 检查是否是历元行
                if len(line.strip()) >= 29 and line[1:3].strip().isdigit():
                    # 检查历元标志
                    flag_str = line[28:29].strip()
                    flag = int(flag_str) if flag_str else 0

                    # 如果历元标志为4，跳过当前行和下一行（COMMENT行）
                    if flag == 4:
                        f.readline()  # 跳过COMMENT行
                        continue

                    # 将文件指针回退一行
                    f.seek(f.tell() - len(line))
                    # 读取观测记录
                    obs_record = self._read_obs_record(f)
                    if not obs_record:
                        break
                    self.obs_data.append(obs_record)

    def _parse_header_line(self, line):
        if "RINEX VERSION / TYPE" in line:
            self.header.version = line[0:20].strip()
            self.header.file_type = line[20:40].strip()
        elif "MARKER NAME" in line:
            self.header.marker_name = line[0:60].strip()
        elif "MARKER NUMBER" in line:
            self.header.marker_number = line[0:60].strip()
        elif "OBSERVER / AGENCY" in line:
            self.header.observer = line[0:20].strip()
            self.header.agency = line[20:60].strip()
        elif "REC # / TYPE / VERS" in line:
            self.header.rec_number = line[0:20].strip()
            self.header.rec_type = line[20:40].strip()
            self.header.rec_version = line[40:60].strip()
        elif "ANT # / TYPE" in line:
            self.header.ant_number = line[0:20].strip()
            self.header.ant_type = line[20:40].strip()
        elif "APPROX POSITION XYZ" in line:
            self.header.approx_pos = [
                float(line[0:14]),
                float(line[14:28]),
                float(line[28:42])
            ]
        elif "ANTENNA: DELTA H/E/N" in line:
            self.header.antenna_delta = [
                float(line[0:14]),
                float(line[14:28]),
                float(line[28:42])
            ]
        elif "# / TYPES OF OBSERV" in line:
            num_obs = int(line[0:6])
            obs_types = line[6:60].strip().split()
            self.header.obs_types = obs_types[:num_obs]
        elif "INTERVAL" in line:
            self.header.interval = float(line[0:10])
        elif "TIME OF FIRST OBS" in line:
            self.header.time_first_obs = [
                int(line[0:6]),  # 年
                int(line[6:12]),  # 月
                int(line[12:18]),  # 日
                int(line[18:24]),  # 时
                int(line[24:30]),  # 分
                float(line[30:43])  # 秒
            ]
        elif "LEAP SECONDS" in line:
            self.header.leap_seconds = int(line[0:6])

    def _read_obs_record(self, f):
        # 读取历元行
        epoch_line = f.readline()
        if not epoch_line:
            return None

        obs = ObsData()

        # 解析历元行
        try:
            # 确保历元行不为空且长度足够
            if len(epoch_line.strip()) < 29:
                return None

            # 解析历元时间
            try:
                year_str = epoch_line[1:3].strip()
                month_str = epoch_line[4:6].strip()
                day_str = epoch_line[7:9].strip()
                hour_str = epoch_line[10:12].strip()
                min_str = epoch_line[13:15].strip()
                sec_str = epoch_line[16:26].strip()

                # 只有当所有时间字段都不为空时才解析
                if all([year_str, month_str, day_str, hour_str, min_str, sec_str]):
                    obs.epoch = [
                        int(year_str),
                        int(month_str),
                        int(day_str),
                        int(hour_str),
                        int(min_str),
                        float(sec_str)
                    ]
                else:
                    return None
            except ValueError:
                print("警告：解析历元时间出错")
                return None

            # 解析历元标志和卫星数
            flag_str = epoch_line[28:29].strip()
            num_sats_str = epoch_line[29:32].strip()

            obs.flag = int(flag_str) if flag_str else 0
            obs.num_sats = int(num_sats_str) if num_sats_str else 0

            # 如果历元标志为4，跳过整个历元
            if obs.flag == 4:
                # 读取下一行（COMMENT行）
                next_line = f.readline()
                return None

            # 读取卫星列表
            sats_str = epoch_line[32:].strip()
            if sats_str:
                obs.sats = [sats_str[i:i + 3].strip() for i in range(0, len(sats_str), 3)]
            else:
                obs.sats = []

            # 读取每颗卫星的观测值
            for sat in obs.sats:
                obs.observations[sat] = {}

                # 计算需要读取的行数
                num_obs_types = len(self.header.obs_types)
                lines_needed = (num_obs_types + 4) // 5  # 每行最多5个观测值
                obs_values = []

                # 读取所有需要的行
                for line_num in range(lines_needed):
                    line = f.readline()
                    if not line:
                        break

                    # 处理该行中的观测值
                    for i in range(0, min(5 * 16, len(line)), 16):
                        val_str = line[i:i + 16]
                        if val_str.strip():
                            try:
                                # 计算当前观测值的类型索引
                                curr_type_idx = line_num * 5 + i // 16
                                if curr_type_idx >= num_obs_types:
                                    continue

                                curr_type = self.header.obs_types[curr_type_idx]

                                # 提取观测值（前14位）
                                value_str = val_str[:14].strip()
                                value = float(value_str) if value_str else 0.0

                                if curr_type in ['L1', 'L2']:
                                    # 对于L1和L2，处理LLI和信号强度
                                    lli_str = val_str[14:15]
                                    strength_str = val_str[15:16]
                                    lli = int(lli_str) if lli_str.strip() else 0
                                    strength = int(strength_str) if strength_str.strip() else 0

                                    obs_values.append({
                                        'value': round(value, 3),
                                        'lli': lli,
                                        'strength': strength
                                    })
                                else:
                                    # 对于其他类型，不处理LLI和信号强度
                                    obs_values.append({
                                        'value': round(value, 3),
                                        'lli': 0,
                                        'strength': 0
                                    })
                            except (ValueError, IndexError) as e:
                                print(f"警告：解析观测值时出错 - {str(e)}")
                                obs_values.append({
                                    'value': 0.0,
                                    'lli': 0,
                                    'strength': 0
                                })
                        else:
                            obs_values.append({
                                'value': 0.0,
                                'lli': 0,
                                'strength': 0
                            })

                # 将观测值与类型对应
                for i, obs_type in enumerate(self.header.obs_types):
                    if i < len(obs_values):
                        obs.observations[sat][obs_type] = obs_values[i]
                    else:
                        obs.observations[sat][obs_type] = {
                            'value': 0.0,
                            'lli': 0,
                            'strength': 0
                        }

        except Exception as e:
            print(f"警告：解析观测记录时出错: {e}")
            return None

        return obs if obs.sats else None

    def save_to_file(self, filename, csv_filename=None):
        """将数据保存到文件
        Args:
            filename: 保存文件头信息的文件名
            csv_filename: 保存观测数据的CSV文件名，如果为None则不保存CSV
        """
        # 保存文件头信息
        with open(filename, 'w', encoding='utf-8') as f:
            f.write("=== 文件头信息 ===\n")
            f.write(f"RINEX版本: {self.header.version}\n")
            f.write(f"文件类型: {self.header.file_type}\n")
            f.write(f"测站名: {self.header.marker_name}\n")
            f.write(f"测站编号: {self.header.marker_number}\n")
            f.write(f"观测者/机构: {self.header.observer} / {self.header.agency}\n")
            f.write(f"接收机信息: {self.header.rec_number} {self.header.rec_type} {self.header.rec_version}\n")
            f.write(f"天线信息: {self.header.ant_number} {self.header.ant_type}\n")
            f.write(f"概略坐标(X/Y/Z): {self.header.approx_pos}\n")
            f.write(f"天线高度(H/E/N): {self.header.antenna_delta}\n")
            f.write(f"观测类型: {self.header.obs_types}\n")
            f.write(f"观测间隔: {self.header.interval}秒\n")
            f.write(f"首次观测时间: {self.header.time_first_obs}\n")
            f.write(f"跳秒数: {self.header.leap_seconds}\n")

        # 保存观测数据到CSV文件
        if csv_filename:
            import csv
            with open(csv_filename, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                # 写入CSV表头
                header = ['历元', '历元标志', '卫星', '观测类型', '观测值', 'LLI', '信号强度']
                writer.writerow(header)

                # 写入观测数据
                for obs in self.obs_data:
                    epoch_str = f"{obs.epoch[0]:02d}/{obs.epoch[1]:02d}/{obs.epoch[2]:02d} " \
                                f"{obs.epoch[3]:02d}:{obs.epoch[4]:02d}:{obs.epoch[5]:06.3f}"

                    for sat in obs.sats:
                        if sat in obs.observations:
                            for obs_type, data in obs.observations[sat].items():
                                row = [
                                    epoch_str,
                                    obs.flag,
                                    sat,
                                    obs_type,
                                    data['value'],
                                    data['lli'] if obs_type in ['L1', 'L2'] else '',
                                    data['strength'] if obs_type in ['L1', 'L2'] else ''
                                ]
                                writer.writerow(row)


# 使用示例
if __name__ == "__main__":
    reader = RinexObsReader()
    reader.read_file(r"C:\Users\xiaocai\Desktop\读取\wuhn1660.12o")  # 替换为你的观测文件名

    # 保存数据到文件
    reader.save_to_file("观测文件头信息.txt", "观测数据.csv")