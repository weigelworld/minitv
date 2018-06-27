# MiniTV
### Alignment viewer using AliTV.

## Setup

```
$ git clone https://github.com/weigelworld/minitv
$ python3 -m virtualenv .minitv-venv
$ source .minitv-venv/bin/activate
$ pip install minitv/
```

`minimap2` must be in your PATH for `minitv` to run.
Either install it globally or download the binary from
the [GitHub repository](https://github.com/lh3/minimap2)
and copy it into `.minitv-venv/bin`.


## Usage

```
minitv --aligner_args "-x asm5 -g 150" --no_tree -a test_data/ACD6/TAIR10.gff3 -a test_data/ACD6/Col-0.gff3 -a test_data/ACD6/KBS-Mac-74.gff3 -a test_data/ACD6/Kn-0.gff3 -a test_data/ACD6/Ty-1.gff3 test_data/ACD6/TAIR10.fa test_data/ACD6/Col-0.fa test_data/ACD6/KBS-Mac-74.fa test_data/ACD6/Kn-0.fa test_data/ACD6/Ty-1.fa > out.alitv.json
```

The resulting `out.alitv.json` file can be loaded into any AliTV instance, such as
[this demo site](https://alitvteam.github.io/AliTV/d3/AliTV.html), by going to the "Import and Export" tab and clicking
the "Load JSON" button.
