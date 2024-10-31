import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

@pytest.fixture
def int_rbs_checker():
    checker = InternalRBSChecker()
    checker.initiate()
    return checker

def test_no_irbs(int_rbs_checker):
    cds = 'TTTGGGAGGAGGTTTATGCCCGGG'
    assert int_rbs_checker.run(cds)[0] == True # because no rbs

def test_yes_irbs(int_rbs_checker):
    cds = 'AGGAGGTTGCGGGTATGCG'
    assert int_rbs_checker.run(cds)[0] == False # has sd seq + M