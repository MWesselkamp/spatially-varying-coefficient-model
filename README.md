# Spatially-varying-coefficient-model

Code is structured as follows.

00 .. : Compute the GAMs
        1: 'ref' gam with thin plate regression splines (trs) but no spatially varying coefficients.
        2: 'basic' gam with spatially varying coeffcients in trs, computed with standardized DF data.
        3: 'ns' gam with spatially varying coeffcients in trs, computed with non-standardized DF data.
        4: 'sl' gam with spatially varying coeffcients in trs, location centered around 0. 
        N: 'neutral' gam with spatially varying coeffcients in trs, computed on a standardized subsample of the DF data. 

01 .. : Examine GAMs
        prepare_data_and_maps: 
        model_coefficient_maps:
        clustering_and_predictions:

02 .. : Niche analysis - Variable-specific effects

03 .. : Niche analysis - GAM Optimium 