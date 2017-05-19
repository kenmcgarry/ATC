# ATC
Anatomical Therapeutic Codes (ATC) are a drug classification system which is extensively used in the field of drug development research. 

Anatomical Therapeutic Codes (ATC) are a drug classification system which is extensively used in the field of drug development research.  There are many drugs and medical compounds that as yet do not have ATC codes, it would be useful to have codes automatically assigned to them by computational methods. Our initial work involved building feedforward multi-layer perceptron models (MLP) but the classification accuracy was poor. To gain insights into the problem we used the Kohonen self-organizing neural network to visualize the relationship between the class labels and the independent variables. The information gained from the learned internal clusters gave a deeper insight into the mapping process. The ability to accurately predict ATC codes was unbalanced due to over and under representation of some ATC classes. Further difficulties arise because many drugs have several, quite different ATC codes because they have many therapeutic uses. We used chemical fingerprint data representing a drug’s chemical structure and chemical activity variables. Evaluation metrics were computed, analysing the predictive performance of various self-organizing  models.

## References
> Chen, F., Jiang, Z.: Prediction of drug's anatomical therapeutic chemical
  (ATC) code by integrating drug-domain data. Journal of Biomedical
  Informatics  58,  80--88 (2015)

> Cheng, X., Zhao, S., Xiao, X., Chou, K.: iATC-mISF: a multi-label classifier
  for predicting the classes of anatomical therapeutic codes. Bioinformatics
  33(3),  341--346 (2016)

> Dunkel, M., Gunther, S., Ahmed, J., Wittig, B.: Superpred: drug classification
  and target prediction. Nucleic Acids Research  36,  W55--W59 (2008)

> Gurulingappa, H., Kolarik, C., Hofmann-Apitius, M., Fluck, J.: Concept-based
  semi-automatic classification of drugs. Journal of Chemical Information
  Modeling  49(8),  1986--1992 (2009)

> Kohonen, T., Oja, E., Simula, O., Visa, A., Kangas, J.: Engineering
  applications of the self-organizing map. Proceedings of the IEEE  84(10),
  1358--1383 (1996)
  
> McGarry, K., Slater, N., Amaning, A.: Identifying candidate drugs for
  repositioning by graph based modeling techniques based on drug side-effects.
  In: The 15th UK Workshop on Computational Intelligence, UKCI-2015. University
  of Exeter, UK (7th-9th September 2015)
