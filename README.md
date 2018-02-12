# SigMat
*A Classification Scheme for Gene Signature Matching*

v1.0.1, updated Feb. 12, 2018.

## Functionality
SigMat is a tool for **matching gene expression signatures to experimental conditions**. For example, given a vector of z-scores of differential expression levels of genes in a sample, SigMat can search through its training drug library to return topK (default: 10) drugs matched to that signature. SigMat can be regarded as a **classification** scheme, since each experimental condition naturally defines a class.

SigMat is able to handle the difficult use case where signatures from a less-studied cell line are matched to signatures from well-studied cell lines. The available but sparse data on a less-studied cell line are utilized by SigMat to "tune" a kernel Support Vector Machine (KSVM) learned on a well-studied cell line. For more algorithmic details, please refer to the **SigMat paper**:

Xiao, J., Blatti, C. and Sinha, S. SigMat: A Classification Scheme for Gene Signature Matching. (Manuscript submitted for review.)

You are welcome to use this tool under the license agreement. Please cite the above paper if SigMat benefits your work.

## Data Preparation
Before runing SigMat, please prepare your data in the format of the sample data under the `data` folder, and replace the sample data files:

- `train/sig_train.csv`: The training data for SigMat. If you do not have your own training data, you may keep the sample training data as it is. The sample training data are a subset of the L1000 GEO dataset ([GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)) processed as described in the SigMat paper. This csv file gives a matrix, whose (*i*,*j*)<sup>th</sup> element represents the differential expression z-score of the *j*<sup>th</sup> gene in the *i*<sup>th</sup> sample. There is no need to provide the gene symbols or sample IDs in this file.
- `train/class_train.csv`: The class of each training sample. In the sample data file, classes are represented by drug IDs. You may leave the sample data file as it is if you do not have your own training data.
- `test/sig_test.csv`: This file can contain either 1) your signature data to be matched to the training library, or 2) a testing data set to evaluate the performance of SigMat. Only in the latter case, you will need to provide `test/class_test.csv`, the ground truth class labels of the testing samples. If so, please make sure that `test/class_test.csv` does not contain any classes not observed in the training data.
- `tune/sig_tune.csv` & `tune/class_tune.csv`: The tuning data for SigMat. SigMat achieves its best signature matching accuracy if provided with tuning signatures from the same cell line as your test signatures. It will work even if your test signatures do not share any common class with your tuning signature! Please make sure that `tune/class_tune.csv` does not contain any classes not observed in the training data.
- `test/class_pred.csv`: This file is not required when you run SigMat. Instead, this is a sample OUTPUT file. Each row gives the top K classes matched to each of your testing signatures. You can specify the value of K when running SigMat.

## Running the Code

To train SigMat and match the testing signatures to the top K classes:

```
cd R_code
Rscript sigmat.R [-k int]
```

To evaluate the performance on the test data and print Success@K to screen:

```
Rscript eval.R
```

The iPython Notebook version of the above two script are also provided. Note that in `sigmat.ipynb` you need to specify the `k` value, which is the `topK` variable.

## Output
`sigmat.R` generates `test/class_pred.csv`. Each row gives the top K classes matched to each of your testing signatures. You can specify the value of K when running SigMat.

`eval.R` evaluates SigMat by comparing `test/class_pred.csv` and `test/class_test.csv` (if provided), and prints the Success@K to screen.
