import pytest
from genedesign.transcript_designer import TranscriptDesigner
@pytest.fixture
def td():
    td = TranscriptDesigner()
    td.initiate()
    return td

def test_phe_freq(td):
    print(td.codonLists['F'])
    assert td.codonLists['F'] == ['TTT']*58 + ['TTC']*42
