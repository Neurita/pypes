import sys
from mock import Mock as MagicMock
from recommonmark.parser import CommonMarkParser


class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
            return Mock()

source_parsers = {
    '.md': CommonMarkParser,
}

source_suffix = ['.rst', '.md']

MOCK_MODULES = ['numpy',
                'scipy',
                'matplotlib', 'matplotlib.pyplot', 'matplotlib.transforms', 'mpl_toolkits.axes_grid1',
                'nipy',
                'hansel',
                'dipy',
                'matplotlib',
                'nipype',
                'kaptan',
                'pydicom',
                'boyle',
                'dcmstack']
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)
