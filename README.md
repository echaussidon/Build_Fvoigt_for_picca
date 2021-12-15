# Compute Fvoigt function for DLA modelisation in the correlation function

**author:** Edmond CHAUSSIDON (LBNL)

**mail:** edmond.chaussidon@cea.fr

**Aim:** This explain how to build a Fvoigt function for picca/fitter2/models/fvoigt_models

## TO CHANGE in build_Fvoigt :

    * DLA catalogue : give the path of a catalogue of DLA

    * QSO catalogue : give the path of a catalogue of quasars

    * Weight lambda : give the path of the weight of each lambda calculated with weight_lambda/build_weight.py

    * output filename : give the path of the output

## How to use weight_lambda/build_weight.py :

    * cf the documentation in build_weight.py

    * the file `data/weight_lambda.txt` is the typical values expected and can be used independently for Saclay/London mocks or DR12.

## OUTPUT :

    * Fvoigt_output.txt file is in the directory output. This input has to be copied in `picca/fitter2/models/fvoigt_models` to be used in *picca*.


## A sucessfull example :

    * (You can ask the file to the author !)

    * Example for Saclay Mocks 4.4.3

    * To change in build_Fvoigt :

            path_dla = 'data/master_DLA_4.4.3.fits'
            path_qso = 'data/zcat_desi_drq_4.4.3.fits'
            path_weight_lambda = 'data/weight_lambda.txt'
            path_output = 'output/Fvoigt_example.txt'

    * Then run in your terminal : (around five minutes for this example)

            `python3 build_Fvoigt.py`
