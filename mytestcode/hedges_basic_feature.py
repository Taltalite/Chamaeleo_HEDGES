import os
import Chamaeleo
from Chamaeleo.utils.pipelines import BasicFeaturePipeline
from Chamaeleo.methods.default import BaseCodingAlgorithm
from Chamaeleo.methods.fixed import Church, Goldman, Grass, Blawat
from Chamaeleo.methods.flowed import DNAFountain, YinYangCode

import sys
sys.path.append('..')
from methods.addition import HEDGES

if __name__ == "__main__":
    root_path = 'D:\workspace\Chamaeleo'
    # root_path = os.path.dirname(Chamaeleo.__file__)

    file_paths = {
        # "Makinami.jpg": os.path.join(root_path,  "mytestcode", "data", "EVA3.png")
        "Mona Lisa.jpg": os.path.join(root_path, "data", "pictures", "Mona Lisa.jpg")
    }
    print('test pic path:')
    for value in file_paths.values():
        print(value)

    coding_schemes = {
        "Base": BaseCodingAlgorithm(),
        "Church et al.": Church(), "Goldman et al.": Goldman(), "Grass et al.": Grass(), "Blawat et al.": Blawat(),
        "DNA Fountain": DNAFountain(redundancy=0.5), "Yin-Yang Code": YinYangCode(),
        "HEDGES": HEDGES()
    }

    needed_indices = [
        True,
        True, True, True, True,
        False, True,
        False       # HEDGES must be false
    ]

    # img = cv2.imread(file_paths.get("Makinami.jpg"))[:,:,(2,1,0)]
    # plt.imshow(img)
    # plt.show()

    pipeline = BasicFeaturePipeline(
        coding_schemes=coding_schemes,
        needed_indices=needed_indices,
        file_paths=file_paths,
        segment_length=100,
        need_logs=True
    )

    pipeline.calculate()
    pipeline.output_records(type="string")

