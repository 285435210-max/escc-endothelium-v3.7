# 发布操作步骤（作者本人执行；本地准备已完成）

本目录已由 `13_refcheck/114_build_code_release_v3_7.py` 生成，只含分析脚本与元数据，
不含原始数据、患者临床表、第三方全文 PDF 与稿件文件。**尚未联网、尚未发布。**

## 1. 账号与绑定（已完成，2026-11-09 实测）
GitHub 账号 `285435210-max`、Zenodo 账号 `285435210@qq.com`、ORCID
0000-0002-0617-3629 均已就绪；Zenodo -> Linked accounts 中 GitHub 与 ORCID
均显示 Disconnect（即已关联）；`/account/settings/github/` 已出现 `Sync now`。
若换浏览器或重新登录，请先确认这三项仍为已连接状态。

## 2. 建公开仓库（约 5 分钟）
- 仓库名建议 `escc-endothelium-v3.7`，**必须 Public**（Zenodo 只能访问公开仓库）。
- 网页方式：New repository -> Public -> 勾选 "Add a README" 不要勾（避免覆盖）；
  创建后用 "Add file -> Upload files" 把本目录内全部文件拖入（含隐藏文件
  `.gitignore` 与 `CITATION.cff`）。
- 或命令行方式（在本目录执行）：
  ```
  git init
  git add -A
  git commit -m "Release v1.0.0: analysis code for manuscript v3.7"
  git branch -M main
  git remote add origin https://github.com/<USER>/escc-endothelium-v3.7.git
  git push -u origin main
  ```
- 上传前自检：仓库内不得出现 `.pdf` / `.docx` / `.xlsx` / `.gz` / `.rds`。

## 3. 绑定 Zenodo 并生成 DOI（约 10 分钟）
1. https://zenodo.org/account/settings/github/ -> 点 `Sync now` -> 新仓库出现后把开关拨到 ON。
2. GitHub 仓库页 -> Releases -> "Create a new release" -> tag 填 `v1.0.0`，
   title 填 `Analysis code for manuscript v3.7` -> Publish release。
3. 等数分钟，Zenodo 自动创建草稿。在 Zenodo -> Uploads 中检查元数据：
   - Upload type: Software
   - Title / Creators（Liu Bin; Ni Zhen; Liang Runbin）/ Affiliation
   - License: MIT
   - Version: 1.0.0
   确认无误后点 Publish。
4. 得到版本 DOI `10.5281/zenodo.XXXXXXX`（指向 v1.0.0）与概念 DOI
   （指向所有版本）。**稿件引用用版本 DOI**。

## 4. 回填稿件
把 DOI 发给助手，由助手一次性改完以下位置并重新生成交付物：
- 英文稿 Methods / Software 段、Data availability statement
- 中文稿对应两段
- 附表 S12_README 与 `S12_analysis_code_inventory`（29 行 Repository status）
- 三份 v3.7 docx、S1-S14 工作簿、独立附表、一致性核验报告

## 注意（发布后的修改规则，2026-11-10 按官方帮助页核对）
- **GitHub 仓库文件**：随时可改，push 新 commit 即可。但公开仓库的旧 commit
  会永久留在历史中，删除不等于未公开。
- **GitHub Release**：可编辑标题与说明，也可删除 release；tag 对应的源码归档
  是固定快照，改代码要让归档更新必须新建 tag 与 Release。若仓库启用了
  immutable releases，发布后连 assets 与 tag 都不能动。
- **Zenodo 元数据**（标题、作者、ORCID、单位、描述、License）：随时可编辑，
  **DOI 不变**。
- **Zenodo 文件**：minor corrections 须在发布后 **30 天内**自助发起
  （Edit published files），编辑草稿须在原发布后 **45 天内**发布，**不改 DOI**；
  超窗仅在版权、个人数据泄露等正当理由下经 support 审批。重大改动应
  New version，生成全新的**版本 DOI**（概念 DOI 不变）。
- **Zenodo 记录删除**：owner 可在发布后 **30 天内**自行删除，之后仅在正当
  理由下申请；已注册 DOI 会保留 tombstone page。
- **稿件引用**：稿件用版本 DOI，指向固定快照。代码后续改版时，需决定是
  回填新 DOI，还是声明稿件对应 v1.0.0。
- **若预计近期还会改代码**：建议先只建 public 仓库并推码，等稿件接收前再发
  v1.0.0 Release 与 Zenodo DOI，避免反复回填稿件。
- 每次新建 GitHub Release 都会生成新的版本 DOI；概念 DOI 保持不变。
- 发布前仍务必再查一遍版权文件。
- 若暂时不做：现有 "available from the corresponding author on reasonable
  request" 仍符合 Frontiers 备用条款，但二审/四审意见中「公开分析代码」一条会保持未响应。
