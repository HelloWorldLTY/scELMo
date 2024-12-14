# <img src="elmo dalle2.png" alt="image_description" width="40" height="40"/>   scELMo: Embeddings from Language Models are Good Learners for Single-cell Data Analysis



# News!

We have uploaded gene embeddings from gpt4-o and drug embeddings from GPT 3.5 in our website, please check them if you wanna have a try!

# Installation

We rely on OpenAI API for query.

```
pip install openai
```

The descriptions and tutorials for OpenAI API can be found in this [link](https://platform.openai.com/).

We reply on these packages for zero-shot learning analysis.

```
pip install scib scib_metrics==0.3.3 pickle mygene scanpy==1.9.3 scikit-learn
```

Installing hnswlib from the original Github profile to avoid potential errors.
```
apt-get install -y python-setuptools python-pip #may not need it for HPC base
git clone https://github.com/nmslib/hnswlib.git
cd hnswlib
pip install .
```
All the packages above are enough for testing tasks absed on zero-shot learning.

We rely on PyTorch for fine-tuning.

```
conda install pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia
conda install lightning -c conda-forge
```

For the perturbation analysis, please install related pacakges based on their website and use the modifeid version provided in the **Perturbation Analysis** folder: [CINEMAOT](https://github.com/vandijklab/CINEMA-OT/tree/main), [CPA](https://github.com/theislab/cpa) and [GEARS](https://github.com/snap-stanford/GEARS/tree/master).

To generate gene embeddings from sequence models (as seq2emb), please refer [seq2cells](https://github.com/GSK-AI/seq2cells) to install related packages. 


For users who cannot access OpenAI API, we provide an alternative solution based on [deepseekv2](https://www.deepseek.com/). Please refer the **Get outputs from LLMs** for more information.

# Tutorials

Please use the example ipynb notebook in each folders as instructions. Evaluations are included in the notebooks. The demo tutorial can be finished in a normal computer within 10 minutes with a prepared environment.

# Datasets

All of the datasets and their download information are included in the Supplementary file 3. A demo dataset for clustering can be found in this [link](https://drive.google.com/file/d/1hHVutJ3tsAhkhTJ-wCNe9OfXubw2m2gN/view?usp=sharing).

# Database for scELMo

We are maintaining a [website](https://sites.google.com/yale.edu/scelmolib) containing embeddings of different information generated by LLM. We are happy to discuss if you have any requests or comments.

# Acknowledgement

We refer the codes from the following packages to implement scELMo. Many thanks to these great developers:

[GenePT](https://github.com/yiqunchen/GenePT), [seq2cells](https://github.com/GSK-AI/seq2cells), [CINEMAOT](https://github.com/vandijklab/CINEMA-OT/tree/main), [CPA](https://github.com/theislab/cpa) and [GEARS](https://github.com/snap-stanford/GEARS/tree/master).

# Open for contribution

We are happy to see if you have more exciting ideas about the extension of scELMo. Feel free to contact us for discussion:

Tianyu Liu (tianyu.liu@yale.edu)

# Citation
```
@article{liu2023scelmo,
  title={scELMo: Embeddings from Language Models are Good Learners for Single-cell Data Analysis},
  author={Liu, Tianyu and Chen, Tianqi and Zheng, Wangjie and Luo, Xiao and Zhao, Hongyu},
  journal={bioRxiv},
  pages={2023--12},
  year={2023},
  publisher={Cold Spring Harbor Laboratory}
}
```

# Related work

- [spEMO](https://github.com/HelloWorldLTY/spEMO)
- [scLAMBDA](https://github.com/gefeiwang/scLAMBDA)