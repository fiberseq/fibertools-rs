```bash
#python3 -m venv .env
source .env/bin/activate
#pip install maturin
#maturin init
#maturin develop
#maturin publish
```


# Remote pip install of another branch
in this case the other branch is called `footprint`.
```bash
pip install -e 'git+https://github.com/fiberseq/fibertools-rs.git@footprint#egg=pyft&subdirectory=py-ft'   
```