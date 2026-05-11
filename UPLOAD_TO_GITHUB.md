# Upload to GitHub

This repository is already configured with the remote:

```text
https://github.com/Echo-nannan/rt-qpcr-protocol-guide.git
```

## Push Future Updates

After editing files locally, run:

```powershell
cd rt-qpcr-protocol-guide
git status
git add .
git commit -m "Update documentation"
git push
```

## First-Time Setup on Another Computer

Clone the repository:

```powershell
git clone https://github.com/Echo-nannan/rt-qpcr-protocol-guide.git
cd rt-qpcr-protocol-guide
```

Install the Python package in editable mode:

```powershell
python -m pip install -e .
```

Run tests:

```powershell
python -m pytest tests -q
```

## Run the Legacy GUI

```powershell
cd legacy\05_Result_Processing\result_processor
python -m pip install -r requirements.txt
python scripts\start_gui.py
```

Generated result folders, logs, virtual environments, executables, and large binary outputs are intentionally ignored by Git.
