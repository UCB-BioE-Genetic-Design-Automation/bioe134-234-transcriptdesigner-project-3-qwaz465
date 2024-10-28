import pytest
from genedesign.transcript_designer import TranscriptDesigner
@pytest.fixture
def td():
    td = TranscriptDesigner()
    td.initiate()
    return td

# has low CAI, what to do here??? better test or get rid of it
def test_get_initial_cds(td):
    init_cds = td.get_initial_cds('WWW')
    assert init_cds == ['TGG', 'TGG', 'TGG']