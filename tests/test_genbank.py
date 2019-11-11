"""
Test suite for genbank.py
"""

from pathlib import Path

import pytest

from clusterblaster import genbank

TEST_DIR = Path(__file__).resolve().parent


def test_get_genbank_paths(tmp_path):
    d = tmp_path / "folder"
    d.mkdir()

    for ext in [".gb", ".gbk", ".genbank", ".fake", ""]:
        _p = d / f"test{ext}"
        _p.write_text("test")

    paths = genbank.get_genbank_paths(d)

    assert paths == [d / "test.genbank", d / "test.gb", d / "test.gbk"]


def test_parse():
    """Test parsing of sample GenBank file in test directory"""
    gbk = TEST_DIR / "sample.gbk"

    with gbk.open() as handle:
        org = genbank.parse(handle)

    assert org == {
        "name": "Saccharomyces cerevisiae",
        "strain": "No_strain",
        "file": str(gbk),
        "scaffolds": [
            {
                "accession": "SCU49845",
                "proteins": [
                    {
                        "index": 0,
                        "id": "AAA98665.1",
                        "start": 1,
                        "end": 206,
                        "strand": "+",
                        "sequence": "SSIYNGISTSGLDLNNGTIADMRQLGIVESYKLKRAVVSSASEAAEVLLRVDNIIRARPRTANRQHM",
                    },
                    {
                        "index": 1,
                        "id": "AAA98666.1",
                        "start": 687,
                        "end": 3158,
                        "strand": "+",
                        "sequence": "MTQLQISLLLTATISLLHLVVATPYEAYPIGKQYPPVARVNESFTFQISNDTYKSSVDKTAQITYNCFDLPSWLSFDSSSRTFSGEPSSDLLSDANTTLYFNVILEGTDSADSTSLNNTYQFVVTNRPSISLSSDFNLLALLKNYGYTNGKNALKLDPNEVFNVTFDRSMFTNEESIVSYYGRSQLYNAPLPNWLFFDSGELKFTGTAPVINSAIAPETSYSFVIIATDIEGFSAVEVEFELVIGAHQLTTSIQNSLIINVTDTGNVSYDLPLNYVYLDDDPISSDKLGSINLLDAPDWVALDNATISGSVPDELLGKNSNPANFSVSIYDTYGDVIYFNFEVVSTTDLFAISSLPNINATRGEWFSYYFLPSQFTDYVNTNVSLEFTNSSQDHDWVKFQSSNLTLAGEVPKNFDKLSLGLKANQGSQSQELYFNIIGMDSKITHSNHSANATSTRSSHHSTSTSSYTSSTYTAKISSTSAAATSSAPAALPAANKTSSHNKKAVAIACGVAIPLGVILVALICFLIFWRRRRENPDDENLPHAISGPDLNNPANKPNQENATPLNNPFDDDASSYDDTSIARRLAALNTLKLDNHSATESDISSVDEKRDSLSGMNTYNDQFQSQSKEELLAKPPVQPPESPFFDPQNRSSSVYMDSEPAVNKSWRYTGNLSPVSDIVRDSYGSQKTVDTEKLFDLEAPEKEKRTSRDVTMSSLDPWNSNISPSPVRKSVTPSPYNVTKHRNRHLQNIQDSQSGKNGITPTTMSTSSSDDFVPVKDGENFCWVHSMEPDRRPSKKRLVDFSNKSNVNVGQVKDIHGRIPEML",
                    },
                    {
                        "index": 2,
                        "id": "AAA98667.1",
                        "start": 3300,
                        "end": 4037,
                        "strand": "-",
                        "sequence": "MNRWVEKWLRVYLKCYINLILFYRNVYPPQSFDYTTYQSFNLPQFVPINRHPALIDYIEELILDVLSKLTHVYRFSICIINKKNDLCIEKYVLDFSELQHVDKDDQIITETEVFDEFRSSLNSLIMHLEKLPKVNDDTITFEAVINAIELELGHKLDRNRRVDSLEEKAEIERDSNWVKCQEDENLPDNNGFQPPKIKLTSLVGSDVGPLIIHQFSEKLISGDDKILNGVYSQYEEGESIFGSLF",
                    },
                ],
            }
        ],
    }
