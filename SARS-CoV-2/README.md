# Minimal set of SARS-CoV-2

Just within two years of COVID-19 pandemic, the number of protein sequence records grows exponentially. For example, in the GISAID EpiCoV, the largest open-access platform for SARS-CoV-2 sequences, the records grew from ~0.02 million to ~125 million (as of Oct 2021) (**Figure 1**). 

![SequenceGISAID_SARS2](https://user-images.githubusercontent.com/51225708/144307768-58242690-8a89-48b8-aa80-8012a6e95027.png)
**Figure 1. The exponential increase in protein sequence records of SARS-CoV-2 in the major public repository, GISAID EpiCoV from January 2020 to October 2021.**

Herein, we demonstrate the application of UNIQmin to SARS-CoV-2 in order to evaluate the efficiency of compression. A total of 56,340,320 (as of July 2021) protein sequence records of SARS-CoV-2 were retrieved from the GISAID EpiCoV database. Consequently, deduplication of the retrieved sequences for each of the 27 SARS-CoV-2 genome encoded proteins resulted in 1,780,901 unique sequences in total (~3.2% of the retrieved dataset). Subsequently, the removal of unique sequences that did not effectively contribute to the peptidome diversity resulted in a minimal set of 273,851 sequences (~0.5% of the retrieved dataset; please refer to the `MinSet_GISAID_SARSCoV2`). This suggests less than ~0.5% of protein sequences of SARS-CoV-2 are alone sufficient to represent the inherent peptidome diversity. 
