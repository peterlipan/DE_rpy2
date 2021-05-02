"""
MIT License

Copyright (c) 2021 Li Pan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import pandas as pd
import numpy as np
import warnings
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri, pandas2ri, Formula
from rpy2.robjects.packages import importr

pandas2ri.activate()
numpy2ri.activate()

# import R libraries
DESeq2 = importr('DESeq2')
edgeR = importr('edgeR')
Limma = importr('limma')
stats = importr('stats')

to_dataframe = robjects.r('function(x) data.frame(x)')


class DE_rpy2:
    """
    Running DESeq2, edgeR, limma through rpy2

    input:
    count_matrix: a pandas dataframe with each column as count
    (float values in FPKM/RPKM are also acceptable as internal rounding will be done)
    , and a id column for gene id
        example:
        id    sampleA    sampleB
        geneA    5.1       1
        geneB    4.2       5
        geneC    1         2
    design_matrix: a pandas dataframe with each column as a condition, and one row for one sample
    Note that the sample name must be the index not a column
                condition
    sampleA1     treated
    sampleA2     treated
    sampleB1    untreated
    sampleB2    untreated

    design_formula: default to be the column name of design matrix, example: "~ condition""
    If it contains multiple conditions, this formula must be customised,
    or the DESeq2 will only consider the first condition.
    gene_column: column name of gene id columns in count_matrix, default = 'id'
    """

    def __init__(self, count_matrix, design_matrix, design_formula=None, gene_column='id'):
        assert gene_column in count_matrix, \
            'column: \'%s\', not found in count matrix' % gene_column
        assert count_matrix.shape[1] - 1 == design_matrix.shape[0], \
            'The number of rows in design matrix must ' \
            'be equal to the number of samples in count matrix'
        assert all(pd.isna(count_matrix)), \
            'Null values are found in count matrix' \
            'Please check it'
        assert len(design_matrix.columns), \
            'Columns names are needed in design matrix'

        if 'float' in count_matrix.drop(gene_column, axis=1):
            warnings.warn('DESeq2 and edgeR only accept integer counts'
                          'The values in count matrix are automatically rounded'
                          'In fact the FPKM/RPKM input is not encouraged by DESeq2 officially')

        # parameters used in DESeq2
        self.count_matrix = pandas2ri.py2ri(count_matrix.drop(gene_column, axis=1).astype('int'))
        self.design_matrix = pandas2ri.py2ri(design_matrix)
        self.gene_ids = count_matrix[gene_column]
        self.gene_column = gene_column
        self.deseq2_result = None
        self.deseq2_label = None
        if design_formula is None:
            condition = design_matrix.columns[0]
            if len(design_matrix.columns) > 1:
                warnings.warn('Multiple conditions are set in design matrix, '
                              'you\'d better customise the design formula.'
                              'Here it only considers the first condition')
            self.design_formula = Formula('~ ' + condition)
        else:
            self.design_formula = Formula(design_formula)

        # parameters used in edgeR
        self.edgeR_group = numpy2ri.py2ri(design_matrix.iloc[:, 0].values)
        self.edgeR_gene_names = numpy2ri.py2ri(count_matrix[gene_column].values)

        self.edgeR_result = None
        self.edgeR_label = None

        # parameters used in limma
        self.limma_result = None
        self.limma_label = None

        self.final_label = None

    def deseq2(self, threshold=0.05, **kwargs):
        """
        Run the standard DESeq2 workflow.
        Get the DESeq2 results as DataFrame.
        Return the label of each gene: 0 for not differentially expressed,
                                       1 for differentially expressed.
        :param threshold: threshold for the adjusted p-value.
                          default = 0.05.
        :param kwargs: parameters of DESeq2 functions.
        See official instructions for details:
        http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
        :return:
        label: pandas.DataFrame format with 2 columns: gene ids and labels
        """
        # Run DESeq2 workflow
        dds = DESeq2.DESeqDataSetFromMatrix(countData=self.count_matrix,
                                            colData=self.design_matrix,
                                            design=self.design_formula)
        dds = DESeq2.DESeq(dds, **kwargs)
        res = DESeq2.results(dds, **kwargs)

        # Store the output matrix as DataFrame
        self.deseq2_result = pandas2ri.ri2py(to_dataframe(res))
        self.deseq2_result[self.gene_column] = self.gene_ids

        # The adjusted p-value in the DESeq2 results
        # may contain NAN
        if any(pd.isna(self.deseq2_result['padj'].values)):
            warnings.warn('There exist NAN in the adjusted p-value'
                          'see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/'
                          'inst/doc/DESeq2.html#why-are-some-p-values-set-to-na')

        # Reject the H0 hypothesis if p-value < threshold
        labels = [int(x) for x in (self.deseq2_result['padj'] < threshold)]
        label = pd.DataFrame({self.gene_column: self.gene_ids, 'label': labels})
        self.deseq2_label = label

        return label

    def edger(self, threshold=0.05):
        """
        Run the standard edgeR workflow.
        Get the edgR results as DataFrame.
        Return the label of each gene:
        0 for not differentially expressed,
        1 for differentially expressed.
        :param threshold: threshold for the p-value.
                          default = 0.05.
        See official instructions for details:
        https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
        :return:
        label: pandas.DataFrame format with 2 columns: gene ids and labels
        """
        # run edgeR workflow
        # Create the DGEList object
        dgList = edgeR.DGEList(counts=self.count_matrix, group=self.edgeR_group, genes=self.edgeR_gene_names)
        # Normalize
        dgList = edgeR.calcNormFactors(dgList, method="TMM")
        # Setting up the model
        robjects.r.assign('edgeR_group', self.edgeR_group)
        designMat = stats.model_matrix(Formula('~ edgeR_group'))
        # Estimating Dispersions
        dgList = edgeR.estimateGLMCommonDisp(dgList, design=designMat)
        dgList = edgeR.estimateGLMTrendedDisp(dgList, design=designMat)
        dgList = edgeR.estimateGLMTagwiseDisp(dgList, design=designMat)
        # Differential Expression
        fit = edgeR.glmQLFit(dgList, designMat)
        test = edgeR.glmQLFTest(fit)
        res = edgeR.topTags(test, n=self.count_matrix.nrow)
        res_df = pandas2ri.ri2py(to_dataframe(res))
        # Sort the result on gene ids
        gene_df = pd.DataFrame({'genes': self.gene_ids})
        self.edgeR_result = pd.merge(gene_df, res_df, how='left')

        # Reject the H0 hypothesis
        labels = [int(x) for x in (self.edgeR_result['PValue'] < threshold)]
        label = pd.DataFrame({self.gene_column: self.gene_ids, 'label': labels})
        self.edgeR_label = label

        return label

    def limma(self, threshold=0.05):
        """
        Run the standard limma workflow.
        Get the limma results as DataFrame.
        Return the label of each gene:
        0 for not differentially expressed,
        1 for differentially expressed.
        :param threshold: threshold for the p-value.
                          default = 0.05.
        See official instructions for details:
        https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
        :return:
        label: pandas.DataFrame format with 2 columns: gene ids and labels
        """
        # Create the DGEList object
        dgList = edgeR.DGEList(counts=self.count_matrix, group=self.edgeR_group, genes=self.edgeR_gene_names)
        # Normalize
        dgList = edgeR.calcNormFactors(dgList, method="TMM")
        # Setting up the model
        robjects.r.assign('edgeR_group', self.edgeR_group)
        designMat = stats.model_matrix(Formula('~ edgeR_group'))
        #  voom
        v = Limma.voom(dgList, designMat)
        # fitting
        fit = Limma.lmFit(v, designMat)
        fit = Limma.eBayes(fit)
        res = Limma.topTable(fit, n=self.count_matrix.nrow)
        res_df = pandas2ri.ri2py(to_dataframe(res))
        # Sort the result on gene ids
        gene_df = pd.DataFrame({'genes': self.gene_ids})
        self.limma_result = pd.merge(gene_df, res_df, how='left')

        # Reject the H0 hypothesis
        labels = [int(x) for x in (self.limma_result['adj.P.Val'] < threshold)]
        label = pd.DataFrame({self.gene_column: self.gene_ids, 'label': labels})
        self.limma_label = label

        return label

    def plot_label_difference(self):
        """
        Plot the Venn diagram of the 3 label output.
        Since we only interest in the differentially expressed genes.
        The number on Venn diagram shows the number of samples labeled as 1.
        Say differentially expressed genes.
        """
        if self.limma_label is None:
            warnings.warn('Seems you haven\'t get limma label'
                          'Automatically running limma...')
            self.limma_label = self.limma()
        if self.deseq2_label is None:
            warnings.warn('Seems you haven\'t get DESeq2 label'
                          'Automatically running DESeq2...')
            self.deseq2_label = self.deseq2()
        if self.edgeR_label is None:
            warnings.warn('Seems you haven\'t get edgeR label'
                          'Automatically running edgeR...')
            self.edgeR_label = self.edger()
        # Import the plot package
        from matplotlib_venn import venn3
        import matplotlib.pyplot as plt

        labels = np.array([self.deseq2_label['label'].values, self.edgeR_label['label'].values,
                           self.limma_label['label'].values]).T
        names = ['DESeq2', 'edgeR', 'limma']
        venn_df = pd.DataFrame(data=labels, columns=names)
        sets = {'000': 0, '001': 0, '010': 0, '011': 0, '100': 0, '101': 0, '110': 0, '111': 0}
        for i in range(venn_df.shape[0]):
            loc = [str(num) for num in venn_df.iloc[i, :]]
            loc = loc[0] + loc[1] + loc[2]
            sets[loc] += 1
        venn3(sets, set_labels=names)
        plt.show()

        return sets

    def get_final_label(self, method='inner'):
        """
        There are 2 methods availabel:
        inner: set those genes as differentially expressed,
        say label 1, if all 3 tools agreed
        vote: set those genes as differentially expressed,
        say label 1, if all 2 out of the 3 tools agreed
        """
        label = None
        menu = ['inner', 'vote']
        assert method in menu, \
            'Please choose the correct method'
        if self.limma_label is None:
            warnings.warn('Seems you haven\'t get limma label'
                          'Automatically running limma...')
            self.limma_label = self.limma()
        if self.deseq2_label is None:
            warnings.warn('Seems you haven\'t get DESeq2 label'
                          'Automatically running DESeq2...')
            self.deseq2_label = self.deseq2()
        if self.edgeR_label is None:
            warnings.warn('Seems you haven\'t get edgeR label'
                          'Automatically running edgeR...')
            self.edgeR_label = self.edger()

        labels = self.deseq2_label['label'].values + self.edgeR_label['label'].values + self.limma_label['label'].values

        if method == 'inner':
            label = [int(x) for x in (labels == 3)]

        if method == 'vote':
            label = [int(x) for x in (labels >= 2)]

        self.final_label = pd.DataFrame({self.gene_column: self.gene_ids, 'label': label})

        return self.final_label
