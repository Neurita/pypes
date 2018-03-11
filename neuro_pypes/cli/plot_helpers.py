import pathlib
from typing import List


def create_imglist_html(img_files: List[pathlib.Path], output_filepath: [str, pathlib.Path]='index.html'):
    """
    """
    if isinstance(output_filepath, str):
        output_filepath = pathlib.Path(output_filepath)

    with output_filepath.open('w') as f:
        f.write('<html>\n')
        f.write('<body bgcolor="#424242">\n')

        for imgf in img_files:
            line = '<p><a href="{imgf}"><img src="{imgf}"/><br/>{imgf}</a></p>'.format(imgf=str(imgf))
            f.write(line + '\n')

        f.write('</body">\n')
        f.write('</html>\n')
