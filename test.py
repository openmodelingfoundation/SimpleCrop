import os
import pytest
import shutil
import subprocess


@pytest.fixture(scope='session', autouse=True)
def initialize(request):
    if not os.path.exists('test/data'):
        shutil.copytree('data', 'test/data')
    root = os.path.dirname(os.path.realpath(__file__))
    subprocess.run([os.path.join(root, 'src/_build/simplecrop')], check=True, cwd=os.path.join(root, 'test'))


def read_soil_lines(fname):
    with open(fname) as f:
        content = f.read()
        lines = content.splitlines()
        return lines[7:]


def read_plant_lines(fname):
    with open(fname) as f:
        content = f.read()
        return content.splitlines()


def test_soil_output():
    soil_orig = read_soil_lines('output/soil.out')
    soil_new = read_soil_lines('test/output/soil.out')
    assert soil_orig == soil_new


def test_plant_output():
    plant_orig = read_plant_lines('output/plant.out')
    plant_new = read_plant_lines('test/output/plant.out')
    assert plant_orig == plant_new