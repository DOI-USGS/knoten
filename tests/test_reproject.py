from unittest import mock
import pytest

from knoten import utils

def test_reproject():
    res = utils.reproject([0,1,0], 10, 10, 'geocent', 'latlon')
    print(res)
    assert res[0] == 0
    assert res[1] == 90
    assert res[2] == -9
