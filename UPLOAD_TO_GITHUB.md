# Upload to GitHub

下面步骤会把当前文件夹上传到你的 GitHub 账号 `Echo-nannan`。

## 方案 A: GitHub 网页上传

1. 打开 <https://github.com/new>
2. Repository name 填写: `rt-qpcr-protocol-guide`
3. Visibility 选择 `Public` 或 `Private`
4. 不要勾选 `Add a README file`，因为本地已经有 README
5. 点击 `Create repository`
6. 进入新仓库后，点击 `uploading an existing file`
7. 把本文件夹 `rt-qpcr-protocol-guide` 里的所有文件拖进去
8. Commit message 填写: `Initial RT-qPCR protocol guide`
9. 点击 `Commit changes`

## 方案 B: Git 命令上传

先在 GitHub 网页创建空仓库：<https://github.com/new>

仓库名建议:

```text
rt-qpcr-protocol-guide
```

然后在 PowerShell 里运行：

```powershell
cd "J:\Learning_Resource\Experiment_Protocal\Basic_Experiment\RT-qpcr\rt-qpcr-protocol-guide"
git init
git add .
git commit -m "Initial RT-qPCR protocol guide"
git branch -M main
git remote add origin https://github.com/Echo-nannan/rt-qpcr-protocol-guide.git
git push -u origin main
```

如果提示没有登录，可以运行：

```powershell
gh auth login
```

按提示选择：

```text
GitHub.com
HTTPS
Login with a web browser
```

登录完成后重新运行：

```powershell
git push -u origin main
```

## 上传前检查

上传前建议确认只包含整理后的干净文件：

```powershell
cd "J:\Learning_Resource\Experiment_Protocal\Basic_Experiment\RT-qpcr\rt-qpcr-protocol-guide"
git status
```

不要把原始 `_archive`、venv、exe、日志和备份目录一起上传。

