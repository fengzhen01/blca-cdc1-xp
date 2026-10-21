# 图题、图注与可直接用于论文的结果文字

## 图6

**中文图题：** GSE293189膀胱癌单细胞转录组的质量控制、无监督聚类与细胞类型注释

**英文图题：** Quality control, unsupervised clustering and cell-type annotation of the GSE293189 bladder cancer single-cell transcriptomic dataset

**图注：** （A）质控后每个细胞检测基因数的分布，红色虚线表示最低基因数阈值。（B）经数据来源批次校正后获得的Leiden无监督聚类UMAP，共识别32个细胞簇。（C）依据典型标志基因注释的主要细胞类型。（D）不同数据来源在校正后UMAP空间中的分布。原始56,192个细胞经质量控制和双细胞去除后保留24,715个细胞。

## 图7

**中文图题：** GSE293189髓系/树突状细胞亚群及交叉呈递评分分布

**英文图题：** Myeloid/dendritic-cell subpopulations and cross-presentation score distribution in GSE293189

**图注：** （A）10,670个髓系/树突状细胞重新降维后形成的Leiden亚群。（B）结合免疫细胞、MHC-II及亚型标志基因识别的高置信cDC1、cDC2、pDC和成熟DC；成熟DC要求LAMP3、CCR7和FSCN1中至少两个基因表达。（C）预先固定的TAP1/TAP2/TAPBP/SEC22B交叉呈递评分在髓系/DC UMAP空间中的分布。肿瘤样本中共识别82个高置信cDC1细胞，来源于14位患者。

## 图S3

**中文图题：** 单细胞主要细胞类型及DC亚群注释的标志基因验证

**英文图题：** Marker-gene validation of major cell-type and dendritic-cell annotations

**图注：** （A）主要细胞类型典型标志基因的表达。点大小表示表达该基因的细胞比例，颜色表示各基因在不同细胞类型间标准化后的平均表达。（B）cDC1、cDC2、pDC、成熟DC及其他髓系细胞的身份和功能标志基因表达。XCR1、CLEC9A和WDFY4支持cDC1注释；CD1C、CLEC10A和FCER1A支持cDC2注释；CLEC4C、IL3RA和GZMB支持pDC注释；LAMP3、CCR7和FSCN1支持成熟DC注释。

## 图8

**中文图题：** 患者配对pseudobulk差异表达及通路富集分析

**英文图题：** Patient-paired pseudobulk differential expression and pathway enrichment analysis

**图注：** （A）8位患者肿瘤样本中cDC1与cDC2的配对pseudobulk差异表达。红色和蓝色分别表示在cDC1中升高和降低且满足FDR<0.10、|log2FC|≥0.5的基因。（B）依据全基因准似然统计量进行的GO/KEGG排序富集分析。蓝色表示相关通路在cDC1中相对降低，即主要偏向cDC2。广义“抗原加工与呈递”通路以MHC-II/溶酶体基因为主，不应等同于cDC1特异的MHC-I交叉呈递程序。

## 图9

**中文图题：** IMvigor210中交叉呈递评分与肿瘤免疫微环境及客观缓解的关系

**英文图题：** Associations of the cross-presentation score with the tumor immune microenvironment and objective response in IMvigor210

**图注：** （A）交叉呈递评分与cDC1身份、T细胞炎症、细胞毒性T细胞、IFN反应、成纤维细胞/TGF-β及免疫抑制性髓系评分的Spearman相关矩阵。（B）缓解者与非缓解者的交叉呈递评分分布。（C）免疫荒漠、免疫排斥和免疫炎症表型间的评分差异。箱体表示四分位距，中线表示中位数，散点表示单个样本。

## 图10

**中文图题：** IMvigor210中交叉呈递评分与疗效及总生存的多变量关联

**英文图题：** Multivariable associations of the cross-presentation score with response and overall survival in IMvigor210

**图注：** 森林图展示交叉呈递评分每升高1个标准差对应的客观缓解优势比（OR）和总生存风险比（HR）及其95%置信区间。模型依次为未校正、校正TMB、校正TMB与cDC1身份，以及校正TMB与T细胞炎症评分。虚线表示无效应值1。

## 可直接用于“结果”部分的凝练文字

对GSE293189的56,192个原始细胞进行质控和样本内双细胞检测后，共保留24,715个细胞。基于高变基因、PCA、邻接图和Leiden算法进行无监督聚类，共获得32个细胞簇，并根据典型标志基因注释为10类主要细胞。主要细胞群在UMAP空间中呈相对清晰的分离，不同数据来源经批次校正后在多数细胞群中发生混合。标志基因点图显示，各主要细胞类型均表达预期的典型标志物，支持聚类和注释结果的可靠性。

进一步对10,670个髓系/DC细胞进行重聚类，并结合免疫细胞、MHC-II和亚型标志基因进行高置信注释。肿瘤样本中共识别82个cDC1细胞，来源于14位患者。cDC1主要表达XCR1、CLEC9A和WDFY4，cDC2主要表达CD1C、CLEC10A和FCER1A，pDC表达CLEC4C、IL3RA和GZMB，成熟DC表达LAMP3、CCR7和FSCN1。固定交叉呈递评分在髓系/DC空间内呈连续分布，提示该功能程序在多个髓系状态之间存在异质性。

为避免细胞级分析造成伪重复，本研究以患者为单位构建pseudobulk计数，并在8位同时含有至少3个cDC1和10个cDC2细胞的患者中开展配对差异分析。cDC1显著高表达XCR1、CLEC9A、WDFY4、DNASE1L3和TAP1，而CLEC10A、IL1B和C1QA相对降低。XCR1、CLEC9A和WDFY4同时用于cDC1注释，其显著差异主要验证了细胞身份的一致性；未参与分类的TAP1仍显著升高，为cDC1中部分交叉呈递程序增强提供了患者级证据。敏感性分析纳入9位患者后差异方向保持一致。

排序通路富集显示，cDC2相对富集脂多糖反应、吞噬、溶酶体、炎症反应及细胞因子–受体相互作用。KEGG“抗原加工与呈递”通路亦整体偏向cDC2，但该通路包含大量MHC-II和溶酶体处理基因，不能直接解释为cDC1交叉呈递减弱。结合TAP1和WDFY4的基因层面结果，更合适的解释是cDC1和cDC2具有不同的抗原处理分工，而当前样本量尚不足以证明完整MHC-I交叉呈递通路在cDC1中整体增强。

在IMvigor210中，固定交叉呈递评分与T细胞炎症、IFN反应和细胞毒性T细胞评分均呈较强正相关，并随肿瘤免疫表型由免疫荒漠向免疫炎症转变而升高。未校正分析显示，高评分与更高客观缓解概率和更长总生存相关；校正TMB后，总生存关联仍显著，而客观缓解关联减弱。进一步校正cDC1身份或T细胞炎症评分后，两类结局的独立关联均不显著。由此表明，交叉呈递评分主要反映cDC1–T细胞炎症轴及免疫活化微环境，目前应视为候选相关标志，而非已经确立的独立预测因子。
