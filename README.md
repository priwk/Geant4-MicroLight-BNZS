# Geant4-MicroLight-BNZS

BN/ZnS 局部微结构 Geant4 模型。

## 运行模式

```text
/cfg/setRunMode StageA_NeutronPatch
/cfg/setRunMode StageB_ReplayAlphaLi
/cfg/setRunMode StageC_OpticalStub
/cfg/setRunMode StageC_OpticalRVE
/cfg/setRunMode StageD_OpticalHomogenization
```

## 源码结构

```text
src/core      通用配置、几何体、物理列表、动作初始化
src/modes     运行模式分发
src/stageA    热中子局部区域输运与 Sigma_eff 汇总
src/stageB    从中子俘获 CSV 记录回放 alpha/Li7
src/stageC    从 ZnS step 源 CSV 记录回放光学光子
src/stageD    随机 RVE 统计均匀化光学参数

include/core
include/modes
include/stageA
include/stageB
include/stageC
include/stageD
```

## Stage A 排布批处理

在一个或多个配比目录下运行所有排布文件：

```bash
STAGEA_EVENTS=100000 bash stageA_batch_placements.sh 1-2 1-3
```

如果不传配比参数，脚本会扫描 `Input/placements` 下的所有目录。
每个配比会输出一个汇总文件：

```text
Output/stageA/<ratio>/neutron_transport_summary.csv
```

Stage A 汇总表会记录 `placement_file`、`placement_basename` 以及排布头中的
`seedBase`，然后为每个排布追加一行 `Sigma_eff` 结果。

Stage C 包含 `StageC_OpticalRVE`，用于基于派生的 ZnS step 源表进行局部光学光子回放。
光学材料属性在几何体中定义，并启用了 `G4OpticalPhysics`。

## Stage D Optical Homogenization

`StageD_OpticalHomogenization` estimates effective optical transport parameters of the random BN/ZnS composite using phase-conditioned statistical re-entry across RVE boundaries. It is not a strict periodic boundary model and does not simulate macroscopic detector thickness. The extracted raw parameters are later ensemble-averaged over multiple placements and calibrated using experimentally measured effective attenuation length.

StageD 使用多套随机 RVE 的同相统计重入来近似无限随机浑浊介质。
它不追求同一个 RVE 的几何周期连续性。
当前主输出口径是 `particle encounter without threshold`：
每次光子进入 BN/ZnS 颗粒时记录入射方向，最终离开该颗粒时记录出射方向，
并据此统计 `mu_a`、`mu_s`、`g1`、`g2`、`mu_s_prime` 和 `P(cos theta)`。
宏观厚度、读出光子、PSF/MTF 由独立 Macro OpticalMC 处理。

第一版实现了：

```text
- RunMode: StageD_OpticalHomogenization
- source_mode: uniform_ZnS
- boundary_mode: same_phase_reentry
- sphere re-entry: same_phase_rho_over_R / same_phase_random
- matrix re-entry: random_matrix
- phase-resolved path length and absorption statistics
- per-photon events CSV
- per-run summary CSV
- tabulated phase function CSV
```

示例宏命令：

```text
/cfg/setRunMode StageD_OpticalHomogenization
/cfg/stageD/setWavelengthNm 450
/cfg/stageD/setSourceMode uniform_ZnS
/cfg/stageD/setBoundaryMode same_phase_reentry
/cfg/stageD/setThetaThresholdDeg 1.0
/cfg/stageD/setMaxReentry 10000
/cfg/stageD/setMaxSteps 100000
/cfg/stageD/setMaxPathLengthUm 1000000
/run/initialize
/run/beamOn 10000
```

输出文件写入：

```text
Output/stageD_optical_homogenization/<ratio>/<placement_id>/lambda_<wavelength>nm/stageD_events.csv
Output/stageD_optical_homogenization/<ratio>/<placement_id>/lambda_<wavelength>nm/stageD_summary.csv
Output/stageD_optical_homogenization/<ratio>/<placement_id>/lambda_<wavelength>nm/phase_function.csv
```

为兼容旧流程，summary 还会额外复制一份：

```text
Output/stageD_optical_homogenization/<ratio>/<placement_id>/lambda_<wavelength>nm/optical_homogenization_summary.csv
```

其中 summary 的主列至少包含：

```text
n_photons
num_absorbed_total
num_absorbed_BN
num_absorbed_ZnS
num_absorbed_Matrix
total_path_length_um
path_length_BN_um
path_length_ZnS_um
path_length_Matrix_um
num_encounter_total
mu_a_count_per_um
mu_a_expected_per_um
mu_s_per_um
g1
g2
mu_s_prime_direct_per_um
mu_s_prime_from_g_per_um
total_reentry
mean_num_reentry
```

StageD 后处理脚本：

```bash
python3 scripts/merge_stageD_optical_params.py --ratio 1-2
python3 scripts/calibrate_optical_params_with_Leff.py \
  --ratio 1-2 \
  --experimental-leff experimental_Leff.csv
```

## Stage C ZnS 源提取

`stageC_make_zns_sources.py` 现在支持两类 Stage B 输入：

- 推荐的 slim 生产输出：
  `*_capture_anchors.csv` + `*_zns_track_steps.csv`
- 兼容调试用的 full 输出：
  `*_alpha_li_steps.csv`

脚本会把 Stage B 的 alpha/Li 轨迹转换成 Stage C 或外部光学 Monte Carlo 可直接使用的 ZnS 源表：

```bash
python3 stageC_make_zns_sources.py
```

脚本会提示输入配比目录，例如 `1-1`、`1-1.5` 或 `1-2`，以及厚度范围，例如 `30-200`。
它会优先扫描：

```text
Output/stageB/<ratio_tag>/*_capture_anchors.csv
Output/stageB/<ratio_tag>/*_zns_track_steps.csv
```

若某个厚度没有 slim 文件对，则回退扫描：

```text
Output/stageB/<ratio_tag>/*_alpha_li_steps.csv
```

并将派生文件写入：

```text
Input/alpha_li_steps/<ratio_tag>/<thickness>_zns_step_sources.csv
Input/alpha_li_steps/<ratio_tag>/<thickness>_event_light_sources.csv
```

只有文件名中的厚度落在指定范围内的文件会被处理，已有厚度会按数值升序处理。

只有 `phase_pre == "ZnS"` 时，才将该 step 的能量沉积计入 ZnS(Ag)。
`phase_post` 会被保留下来，用于边界或跨相位诊断，但不会用于统计 ZnS 发光。
后续 Stage C 跟踪默认使用的光学源模型，是沿每个保留的 pre-to-post step 线段做均匀采样；
派生 step 表中的中点坐标只是调试字段。

在 slim 模式下：

- `capture_anchors.csv` 的每一行代表一条 alpha/Li 轨迹样本
- `physical_event_uid = eventID_record_index`
- `source_event_uid = eventID_record_index_replayIndex`
- `trajectory_weight = 1 / alphali_replay_count`

`event_light_sources.csv` 基于所有 Stage B trajectory anchor 构建，再与 ZnS step 汇总结果做左连接。
即使某条 trajectory 没有 ZnS 能量沉积，也会保留下来，并令 `n_photon0 = 0`，
从而避免平均发光量偏向非零发光 trajectory。

`phase_pre == "ZnS"` 的 step 才会被视为 ZnS 潜在发光轨迹；
`phase_post` 只作为边界/跨相位诊断元数据保留。

可选的闪烁与 Birks 参数：

```bash
python3 stageC_make_zns_sources.py \
  --yield-zns 55000 \
  --birks-kb 0
```

如果要进行 Stage C 光学输运，可按排布拆分派生源表，使每次 Geant4 运行只使用一个 RVE 几何体：

```bash
python3 stageC_make_zns_sources.py \
  --ratio-tag 1-1.5 \
  --thickness-range 500-500 \
  --split-by-placement
```

注意：当前 `--split-by-placement` 仍用于 legacy full 输入
`*_alpha_li_steps.csv` 的拆分流程。
slim 输出当前主要面向后处理或外部光学 MC，
以及不拆分排布的 `stageC_make_zns_sources.py` 批量转换流程。

拆分模式会写出类似下面的文件：

```text
Input/alpha_li_steps/<ratio_tag>/by_placement/<thickness>/<thickness>_<placement_stem>_zns_step_sources.csv
Input/alpha_li_steps/<ratio_tag>/by_placement/<thickness>/<thickness>_<placement_stem>_event_light_sources.csv
```

## Stage C Optical RVE 运行模式

`StageC_OpticalRVE` 会复用一个排布 RVE，并读取一个派生的单排布
`<thickness>_<placement_stem>_zns_step_sources.csv`。
它不会重新生成 alpha/Li 轨迹。
每个 Geant4 事件都会从一个保留的 ZnS step 中采样并发射一个 450 nm 光学光子。
默认源模型会沿原始的 pre-to-post step 线段均匀采样；`midpoint` 只作为调试捷径提供。

示例宏命令：

```text
/cfg/setRunMode StageC_OpticalRVE
/cfg/setOpticalSourceCsv ../Input/alpha_li_steps/1-2/by_placement/30/30_placement_f_0.64_0001_zns_step_sources.csv
/cfg/setSourceSampling uniformAlongStep
/cfg/setSamplePhotonsPerStep 10
/cfg/setOpticalParams 1.5 1000000 2.1 10 2.36 50
/run/initialize
/run/beamOn 2000000
```

RVE 排布会从源 CSV 第一行数据中的 `placement_file` 列自动推断。
除非 `/cfg/setPlacementFile` 与该排布文件完全一致，否则不要手动设置它；
现在 Stage C 会拒绝不匹配的情况。
如果源 CSV 后续行包含不同的 `placement_file`，Stage C 也会拒绝该文件，
因为一次运行不能混用多个 RVE 排布。

对每个源 step，所有采样光子都携带：

```text
photon_weight = n_photon_step / N_sample_photons_per_step
```

输出文件写入：

```text
Output/stageC/<ratio_tag>/<source_stem>_optical_photons.csv
Output/stageC/<ratio_tag>/<source_stem>_optical_summary.csv
Output/stageC/<ratio_tag>/<source_stem>_local_optical_kernel_events.csv
Output/stageC/<ratio_tag>/<source_stem>_local_optical_kernel_steps.csv
Output/stageC/<ratio_tag>/<source_stem>_local_optical_exit_photons.csv
```

默认情况下，为了避免逐光子 CSV 导致输出体积过大，
`*_optical_photons.csv` 不再写出。
如需调试逐光子结局，可显式开启：

```text
/cfg/setWriteStageCPhotonCsv true
```

或设置环境变量：

```bash
BNZS_WRITE_STAGEC_PHOTON_CSV=1
```

`*_local_optical_kernel_events.csv` 会按原始 Stage B 的 `eventID` 聚合采样光子权重，
报告前表面、后表面、侧面逃逸权重，吸收/损失权重，逃逸效率，平均逃逸角，
以及平均路径长度。`*_local_optical_exit_photons.csv` 只保留逃逸光子，
用于出口位置和角度诊断。

## Stage C 宏观耦合

将局部 RVE 光学输出耦合到一个宏观衰减长度：

```bash
python3 stageC_macro_coupling.py \
  --ratio-tag 1-1.5 \
  --l-att-um 1000 \
  --macro-model depth-only
```

如果要使用拟合得到的 UV-Vis 或实验衰减输入，可使用 `--attenuation-csv`：

```text
ratio_label,bn_wt,zns_wt,L_att_um,mu_eff_per_um,source,fit_method
1-1.5,1,1.5,1000,0.001,placeholder,manual
```

脚本会写出：

```text
Output/stageC/<ratio_tag>/thickness_light_yield_curve.csv
```

如果要额外导出可用于后续 PSF / MTF / 分辨率分析的“读出面光斑点表”，
可加上：

```bash
python3 stageC_macro_coupling.py \
  --ratio-tag 1-1.5 \
  --l-att-um 1000 \
  --macro-model angle-resolved \
  --spot-output
```

这会额外写出：

```text
Output/stageC/<ratio_tag>/readout_spot_table_<macro_model>.csv
```

该表会为每个参与宏观传播的逃逸光子保留：

- 事件锚点：`source_event_uid`, `eventID`, `placement_file`
- 宏观俘获位置：`capture_x_um`, `capture_y_um`, `capture_depth_um`
- 局部锚点：`local_capture_x_um`, `local_capture_y_um`, `local_capture_z_um`
- 局部出口态：`exit_surface`, `exit_x_um`, `exit_y_um`, `exit_z_um`, `exit_dir_x/y/z`
- 全局出口位置：`exit_global_x_um`, `exit_global_y_um`, `exit_global_depth_um`
- 读出面位置与权重：`readout_x_um`, `readout_y_um`, `macro_path_um`, `readout_weight`

`depth-only` 对前端读出使用 `s = capture_depth_um`，
对后端读出使用 `s = thickness_um - capture_depth_um`。
`angle-resolved` 使用 `s / max(abs(exit_dir_z), epsilon)`。
当前宏观耦合会继续利用 `*_local_optical_exit_photons.csv` 中的侧向逃逸光子：
先结合 Stage B 事件锚点恢复全局出口位置，再沿出口方向投影到前/后读出面，
并乘以宏观衰减。输出中的
`mean_macro_side_to_front_light_per_capture` 与
`mean_macro_side_to_back_light_per_capture`
会单独报告这部分贡献。

默认情况下，`stageC_macro_coupling.py` 会从外部宏观模型结果读取厚度级捕获统计：

```text
Input/stageA/<ratio_tag>/neutron_transport_summary/neutron_transport_summary.csv
```

优先使用该表中的 `n_incident` 与 `n_absorb` 计算 `P_capture` 和
`*_per_incident_neutron` 列；如果缺失，再回退到 `sigma_eff_per_um`。

如果要在不手工编辑宏文件的情况下运行多个拆分排布，先生成一批 Stage C 宏：

```bash
python3 stageC_run_batch.py \
  --ratio-tag 1-1.5 \
  --thickness 500 \
  --max-placements 5 \
  --samples-per-step 1
```

这些宏会写入：

```text
Output/stageC_macros/<ratio_tag>/<thickness>/
```

加上 `--run` 可对每个生成的宏执行 `build/Geant4-MicroLight-BNZS`。
加上 `--run-macro-coupling --l-att-um <value>`，
可在所选排布运行结束后继续执行宏观耦合步骤。

如果要按完整厚度批量运行，可使用：

```bash
python3 stageC_run_thickness_batch.py \
  --ratio-tag 1-1.5 \
  --thicknesses 30,50,100,200,500,700,1000 \
  --samples-per-step 1 \
  --max-placements 0 \
  --l-att-um 1000 \
  --bootstrap 500
```

建议先加上 `--dry-run` 预览命令，再正式运行 Geant4。

在新增 Stage C 源文件后，构建前请重新运行 CMake，
这样生成的构建系统才会包含 `src/stageC/*.cc` 文件。

## Stage B 平衡 alpha/Li 回放批处理

Stage B 现在默认优先从以下位置读取俘获记录：

```text
Input/stageA/<ratio>/neutron_capture_positions/*_neutron_capture_positions.csv
```

向后兼容的旧路径仍保留为回退：

```text
Input/neutron_capture_positions/<ratio>/*_neutron_capture_positions.csv
```

使用平衡的排布轮转方式运行所有配比目录：

```bash
bash batch_run.sh
```

运行指定的配比目录：

```bash
bash batch_run.sh 1-2 1-3
```

该批处理会先打乱每个兼容的俘获 CSV，然后将其中记录平均拆分到以下目录中的所有排布文件：

```text
Input/placements/<ratio>/*.csv
```

每次 `Geant4-MicroLight-BNZS` 调用只使用一个固定排布和一个临时俘获分块，因此 Geant4 单次运行期间几何体不会变化。

这里有两种不同的“重复”语义：

1. `STAGEB_REPLAY_MULTIPLIER`
   让同一批 capture records 以不同的排布轮转方式重复批处理分配。
2. `BNZS_ALPHALI_REPLAY_PER_CAPTURE`
   对每个物理中子俘获点，在同一个 `local_capture_x/y/z` 上重复发射 `N` 次 alpha/Li，
   每次重新随机反应分支与发射方向。

如果要在额外的旋转排布分配上重复回放同一批已打乱的俘获记录：

```bash
STAGEB_REPLAY_MULTIPLIER=3 bash batch_run.sh 1-2
```

如果要对每个俘获点做两次条件 alpha/Li trajectory 采样：

```bash
BNZS_ALPHALI_REPLAY_PER_CAPTURE=2 bash batch_run.sh 1-2
```

推荐的生产模式输出：

```bash
BNZS_STAGEB_OUTPUT_MODE=slim bash batch_run.sh 1-2
```

slim 模式会按配比和厚度写入：

```text
Output/stageB/<ratio>/<thickness>_capture_anchors.csv
Output/stageB/<ratio>/<thickness>_zns_track_steps.csv
Output/stageB/<ratio>/<thickness>_boundary_stop_summary.csv
```

其中：

- `physical_event_uid = eventID_record_index`
- `source_event_uid = eventID_record_index_replayIndex`
- `trajectory_weight = 1 / alphali_replay_count`

`capture_anchors.csv` 保留所有 replay trajectory，
即使某条 trajectory 没有任何 ZnS step；
`zns_track_steps.csv` 只保留 `phase_pre == "ZnS"` 且 `step_len_um > 0` 的 alpha/Li 轨迹段。

legacy 调试模式 `BNZS_STAGEB_OUTPUT_MODE=full` 仍会输出：

```text
Output/stageB/<ratio>/<thickness>_alpha_li_steps.csv
Output/stageB/1-2/400_alpha_li_steps.csv
```

平衡轮转分块输出的文件名会包含源厚度、回放轮次和排布标签，例如：

```text
Output/stageB/1-2/400_m01_p0001_placement_f_0.64_0001_alpha_li_steps.csv
```

当某一厚度下的所有排布分块都完成后，它们会被合并为正式分析文件，
并且这些分块文件默认会被删除：

```text
Output/stageB/<ratio>/<thickness>_capture_anchors.csv
Output/stageB/<ratio>/<thickness>_zns_track_steps.csv
Output/stageB/<ratio>/<thickness>_boundary_stop_summary.csv
```

如需保留分块输出以便调试：

```bash
python3 stageB_balanced_cycle.py --keep-part-outputs 1-2
```

如需在不重新运行 Geant4 的情况下合并已生成的分块输出：

```bash
python3 stageB_balanced_cycle.py --merge-only 1-2
```

Stage B 会记录真实的局部 alpha/Li 生成点：
`local_capture_x_um`、`local_capture_y_um` 和 `local_capture_z_um`，
以及该次运行使用的 `placement_file`。
表面俘获会保留从宏观前后表面距离到局部 RVE z 切片的映射。
体内俘获会将 `depth_um` 作为元数据保留，并在带有 `6.8 um` XY 安全边距和 `7.4 um` z 安全边距的块体 BN 体积中采样局部俘获点。

日志文件会包含配比和排布基础名：

```text
logs/stageB/balanced/<ratio>/<thickness>/mXX/pYYYY_<placement>.log
```

## 从零重跑一套 Stage B 到 Stage C 数据

如果你想基于当前代码重新生成一套完整的 `Stage B` 与 `Stage C` 数据，
推荐按下面的顺序执行：

1. 先重跑 `Stage B`，推荐生成新的 slim 输出
2. 再用新的 `Stage B` 输出重建 `Stage C` 源表
3. 再运行 `Stage C Optical RVE`
4. 最后做 `Stage C` 宏观耦合曲线

这样可以避免继续沿用旧的 `alpha_li_steps.csv`、旧的 slim/full 混合结果、旧的拆分光学源表，
或者把新旧结果混在一起。

### 0. 建议先清理或备份旧结果

如果你要对某个配比做“干净重跑”，建议先手动备份或删除该配比对应的旧结果目录。
例如，对 `1-2` 配比，常见的相关目录有：

```text
Output/stageB/1-2
Input/alpha_li_steps/1-2
Output/stageC/1-2
Output/stageC_macros/1-2
logs/stageB/balanced/1-2
```

如果不清理，脚本通常仍然能运行，但目录中可能会同时存在旧文件和新文件，
后处理时更容易混淆。

### 1. 重跑 Stage B

最常用的入口脚本是：

```bash
cd /home/jagora/g4work/B2
bash batch_run.sh 1-2
```

如果不带参数，脚本会扫描 `Input/placements/` 下的所有配比目录：

```bash
bash batch_run.sh
```

如果想一次跑多个配比：

```bash
bash batch_run.sh 1-2 1-3
```

这个脚本会自动：

1. 进入 `build/`
2. 运行 `cmake ..`
3. 运行 `make -j$(nproc)`
4. 调用 `stageB_balanced_cycle.py`

`stageB_balanced_cycle.py` 会把每个厚度的俘获记录打乱，再平均分配到同一配比目录下的多个排布文件，
然后每次调用一次 `Geant4-MicroLight-BNZS`，只使用一个固定排布和一个临时俘获分块。

常用的 `Stage B` 额外参数：

```bash
STAGEB_REPLAY_MULTIPLIER=3 bash batch_run.sh 1-2
STAGEB_SHUFFLE_SEED=20260427 bash batch_run.sh 1-2
```

- `STAGEB_REPLAY_MULTIPLIER`
  会让同一批已打乱的俘获记录，以不同的排布轮转方式重复回放。
- `STAGEB_SHUFFLE_SEED`
  控制打乱与排布轮转的随机种子，便于复现。

`Stage B` 主要输入目录：

```text
Input/stageA/<ratio>/neutron_capture_positions/*_neutron_capture_positions.csv
Input/placements/<ratio>/*.csv
```

`Stage B` 主要输出目录：

```text
Output/stageB/<ratio>/<thickness>_capture_anchors.csv
Output/stageB/<ratio>/<thickness>_zns_track_steps.csv
Output/stageB/<ratio>/<thickness>_boundary_stop_summary.csv
logs/stageB/balanced/<ratio>/<thickness>/mXX/pYYYY_<placement>.log
```

如果你想直接调用底层脚本，也可以：

```bash
python3 stageB_balanced_cycle.py 1-2
```

常用选项：

```bash
python3 stageB_balanced_cycle.py \
  --replay-multiplier 2 \
  --seed 20260427 \
  --min-thickness-um 30 \
  --keep-part-outputs \
  1-2
```

其中：

- `--replay-multiplier`：等价于环境变量 `STAGEB_REPLAY_MULTIPLIER`
- `--seed`：等价于环境变量 `STAGEB_SHUFFLE_SEED`
- `--min-thickness-um`：跳过小于该值的厚度文件
- `--keep-part-outputs`：保留分块输出，不只保留合并后的正式文件
- `--merge-only`：只合并已有分块，不重新运行 Geant4
- `--dry-run`：打印将执行的运行配置，不真正调用 `Geant4-MicroLight-BNZS`

常用环境变量：

- `BNZS_STAGEB_OUTPUT_MODE=slim`
  推荐生产输出模式
- `BNZS_ALPHALI_REPLAY_PER_CAPTURE=2`
  对同一物理俘获点做两次 alpha/Li 条件 trajectory 采样
- `STAGEB_THICKNESSES=30,40,...`
  只处理指定厚度的 capture CSV

推荐的生产运行示例：

```bash
STAGEB_THICKNESSES=30,40,50,70,100,125,150,175,200,225,250,275,300,325,350,375,400,500,700,1000 \
BNZS_STAGEB_OUTPUT_MODE=slim \
bash batch_run.sh 1-2
```

如果希望对每个物理中子俘获点做两次 alpha/Li 条件 trajectory 采样：

```bash
STAGEB_THICKNESSES=30,40,50,70,100,125,150,175,200,225,250,275,300,325,350,375,400,500,700,1000 \
BNZS_STAGEB_OUTPUT_MODE=slim \
BNZS_ALPHALI_REPLAY_PER_CAPTURE=2 \
bash batch_run.sh 1-2
```

批处理脚本会为每个临时分块自动生成匹配的 `/run/beamOn` 数量：

```text
beamOn = capture records in chunk × BNZS_ALPHALI_REPLAY_PER_CAPTURE
```

### 2. 用新的 Stage B 输出重建 Stage C 源表

`Stage C` 的局部光学输入不是直接读光学光子，
而是先经过 `stageC_make_zns_sources.py` 把 Stage B 的 trajectory 输出转成两类派生源表：

```text
Input/alpha_li_steps/<ratio>/<thickness>_zns_step_sources.csv
Input/alpha_li_steps/<ratio>/<thickness>_event_light_sources.csv
```

如果要给 `Stage C Optical RVE` 使用，通常还要按排布拆分：

```text
Input/alpha_li_steps/<ratio>/by_placement/<thickness>/<thickness>_<placement>_zns_step_sources.csv
Input/alpha_li_steps/<ratio>/by_placement/<thickness>/<thickness>_<placement>_event_light_sources.csv
```

对单个配比、单段厚度范围，推荐这样做：

```bash
python3 stageC_make_zns_sources.py \
  --ratio-tag 1-2 \
  --thickness-range 30-1000 \
  --split-by-placement
```

常用参数：

```bash
python3 stageC_make_zns_sources.py \
  --ratio-tag 1-2 \
  --thickness-range 30-1000 \
  --split-by-placement \
  --yield-zns 55000 \
  --birks-kb 0
```

- `--yield-zns`
  控制 ZnS(Ag) 的单位能量发光产额，单位是 `photons / MeV`
- `--birks-kb`
  控制 Birks 淬灭；`0` 表示关闭淬灭
- `--split-by-placement`
  为每个 `placement_file` 单独输出一对源表，这是 `StageC_OpticalRVE` 的推荐输入形式

如果你已经确认 `Stage B` 输出没有变化，只是想重复跑后续光学部分，
才适合复用旧源表。
如果 `Stage B` 已经重跑，通常不要加 `--skip-source`，而应重新生成一遍源表。

### 3. 运行 Stage C Optical RVE

最方便的批量入口是：

```bash
python3 stageC_run_thickness_batch.py \
  --ratio-tag 1-2 \
  --thicknesses 30,40,50,60,70,80,90,100,125,150,175,200,225,250,275,300,325,350,375,400,500,700,1000 \
  --samples-per-step 1 \
  --max-placements 0 \
  --l-att-um 1000 \
  --optical-params "1.5 1000000 2.1 10 2.36 50"
```

这条命令会依次做：

1. 对每个厚度调用 `stageC_make_zns_sources.py --split-by-placement`
2. 对每个厚度调用 `stageC_run_batch.py --run`
3. 在全部厚度完成后调用 `stageC_macro_coupling.py`

其中：

- `--ratio-tag`
  指定配比目录，例如 `1-2`
- `--thicknesses`
  指定要跑的厚度列表，逗号或空格分隔
- `--samples-per-step`
  每个保留的 ZnS step 采样多少个光学光子事件
- `--max-placements`
  每个厚度最多跑多少个排布；`0` 表示全部排布
- `--l-att-um`
  宏观耦合时使用的衰减长度，单位 `um`
- `--optical-params`
  六个光学参数，顺序是：
  `matrix_n matrix_abs_um bn_n bn_abs_um zns_n zns_abs_um`

如果想先看命令，不实际运行：

```bash
python3 stageC_run_thickness_batch.py \
  --ratio-tag 1-2 \
  --thicknesses 30,40,50,60,70,80,90,100,125,150,175,200,225,250,275,300,325,350,375,400,500,700,1000 \
  --samples-per-step 1 \
  --max-placements 0 \
  --l-att-um 1000 \
  --optical-params "1.5 1000000 2.1 10 2.36 50" \
  --dry-run
```

如果你已经重新生成好了 `Stage C` 源表，只想跳过第 2 步源表生成：

```bash
python3 stageC_run_thickness_batch.py \
  --ratio-tag 1-2 \
  --thicknesses 30,50,100,200,500,700,1000 \
  --samples-per-step 1 \
  --max-placements 0 \
  --l-att-um 1000 \
  --optical-params "1.5 1000000 2.1 10 2.36 50" \
  --skip-source
```

但要注意：只有在 `Stage B` 输出没变化时，`--skip-source` 才安全。

### 4. 手动拆开跑 Stage C 的三步

如果你不想用总控脚本，而是想把 `Stage C` 分成“生成源表、运行光学、宏观耦合”三步手动执行：

第一步，生成按排布拆分的源表：

```bash
python3 stageC_make_zns_sources.py \
  --ratio-tag 1-2 \
  --thickness-range 500-500 \
  --split-by-placement
```

第二步，只跑某一个厚度的光学回放：

```bash
python3 stageC_run_batch.py \
  --ratio-tag 1-2 \
  --thickness 500 \
  --samples-per-step 1 \
  --max-placements 0 \
  --optical-params "1.5 1000000 2.1 10 2.36 50" \
  --run
```

第三步，对已生成的局部光学 kernel 做宏观耦合：

```bash
python3 stageC_macro_coupling.py \
  --ratio-tag 1-2 \
  --l-att-um 1000 \
  --macro-model depth-only
```

如果你还想输出 angle-resolved 曲线：

```bash
python3 stageC_macro_coupling.py \
  --ratio-tag 1-2 \
  --l-att-um 1000 \
  --macro-model angle-resolved
```

### 5. Stage C 运行后的主要输出

局部光学回放输出位于：

```text
Output/stageC/<ratio>/<source_stem>_optical_photons.csv
Output/stageC/<ratio>/<source_stem>_optical_summary.csv
Output/stageC/<ratio>/<source_stem>_local_optical_kernel_events.csv
Output/stageC/<ratio>/<source_stem>_local_optical_kernel_steps.csv
Output/stageC/<ratio>/<source_stem>_local_optical_exit_photons.csv
```

其中：

- `*_optical_photons.csv`
  默认关闭；仅在 `/cfg/setWriteStageCPhotonCsv true` 或
  `BNZS_WRITE_STAGEC_PHOTON_CSV=1` 时写出，用于记录每个采样光子的结局、
  终点位置和权重
- `*_local_optical_exit_photons.csv`
  只保留成功从前、后、侧面逃逸的光子
- `*_local_optical_kernel_events.csv`
  按原始 `Stage B eventID` 聚合局部 kernel
- `*_local_optical_kernel_steps.csv`
  按原始源 step 聚合局部 kernel
- `*_optical_summary.csv`
  给出该源文件整体的前/后/侧逃逸权重比例

宏观耦合后的厚度曲线输出位于：

```text
Output/stageC/<ratio>/thickness_light_yield_curve.csv
Output/stageC/<ratio>/thickness_light_yield_curve_angle.csv
```

其中第二个文件只有在启用了 `angle-resolved` 曲线时才会生成。

### 6. 常见建议

- 建议优先使用 `bash batch_run.sh <ratio>` 重跑 `Stage B`，
  因为它会自动编译并调用平衡回放脚本。
- 只要 `Stage B` 重跑过，就建议同时重建 `Input/alpha_li_steps/<ratio>`，
  不要继续沿用旧的拆分源表。
- `Stage C` 总控脚本第一次正式跑前，建议先加一次 `--dry-run` 检查命令和厚度列表。
- 如果只是测试流程，可以把 `--max-placements` 设成较小值，例如 `1` 或 `5`，
  先确认输出结构和日志都正常，再全量运行。
