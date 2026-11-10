# BLCA-cDC1-XP：膀胱癌cDC1交叉呈递分析

**BLCA-cDC1-XP** 是 **Bladder Cancer cDC1 Cross-Presentation** 的缩写。本仓库公开膀胱癌cDC1身份、抗原交叉呈递程序及免疫检查点治疗结局分析的可复现代码、必要派生数据、统计结果和论文级图件。

> 当前状态：论文准备阶段的分析仓库，仅用于科研复现，不属于临床诊断或治疗决策工具。

## 研究问题

本项目主要回答三个问题：

1. 能否在膀胱癌单细胞数据中稳健识别cDC1及其他树突状细胞状态；
2. 患者层面的cDC1与cDC2转录程序存在哪些差异；
3. 预先锁定的交叉呈递评分是否与炎症型肿瘤微环境、免疫治疗缓解及总生存相关。

固定交叉呈递评分由 **TAP1、TAP2、TAPBP和SEC22B** 构成。XCR1、CLEC9A和WDFY4等cDC1身份标志单独分析，以减少定义循环。

## 核心结果

- GSE293189标准流程保留24,715个细胞，获得32个Leiden簇；
- 髓系/DC细胞共10,670个，在14位患者的肿瘤样本中识别82个高置信cDC1细胞；
- 8位患者的配对pseudobulk主分析显示，cDC1高表达XCR1、CLEC9A、WDFY4、DNASE1L3和TAP1；
- 广义抗原加工与呈递通路并非整体偏向cDC1，MHC-II/溶酶体程序更多反映cDC2功能，因此结果支持部分交叉呈递基因增强，而非完整通路全面激活；
- IMvigor210中交叉呈递评分与免疫炎症状态、客观缓解和总生存相关，但校正TMB、cDC1身份或T细胞炎症状态后独立关联减弱。

目前最稳妥的定位是：该评分反映**cDC1–T细胞炎症轴**，尚不能作为已确立的独立临床预测标志物。

## 使用方法

安装环境：

```bash
conda env create -f environment.yml
conda activate blca-cdc1-xp
```

原始数据获取和目录结构见 [docs/DATA.md](docs/DATA.md)，完整复现顺序见 [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md)，结果边界见 [docs/RESULTS_INTERPRETATION.md](docs/RESULTS_INTERPRETATION.md)。

## 数据说明

GitHub仓库不直接上传以下内容：原始单细胞表达矩阵、完整H5AD对象、IMvigor210原始R包数据、第三方论文PDF及本地Python/R依赖目录。仓库提供公开数据获取方式、必要的患者级派生数据和统计结果表，从而兼顾可复现性、仓库体积和数据许可要求。

## 许可与声明

代码采用MIT许可证。第三方数据和论文仍受原始许可约束。本项目仅供科研使用，不得用于临床诊断、用药选择或患者管理。
