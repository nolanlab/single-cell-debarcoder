%% this is a sample script that debarcodes an fcs file using the scd.m class without using the GUI

obj = scd('sample_barcode_key.csv');
obj = obj.load_fcs_files('sample_barcoded_file.fcs');

obj = obj.find_bc_cols_by_mass;
obj = obj.load_bcs;

%% normalize each population based on preliminary debarcoding

obj=obj.compute_debarcoding('bcs');
obj=obj.normalize_by_pop;
obj=obj.compute_debarcoding;

%% replace above with this section to optionally normalize by adjusting cofactors for each barcoding channel

% obj = obj.normalize_bcs('bcs');
% obj = obj.compute_debarcoding;
% obj=obj.calculate_cofactors; 
% obj=obj.recofactor; %calculates cofactored_bcs from bcs and cofactors
% obj=obj.normalize_bcs('cofactored_bcs'); %calculates normbcs from cofactored_bcs
% obj=obj.compute_debarcoding;

%%

obj.sep_cutoff = 0.35; %separation threshold

obj = obj.compute_mahal; %

obj.mahal_cutoff_val = 30; %mahalanobis threshold

%write fcs files for each barcode sample and an unassigned file
obj.write_bc_fcs_files('tmp','test')

%compute vector of assigned barcode labels
bclabels=obj.bcind;
bclabels(obj.deltas<=obj.sep_cutoff | obj.mahal>=obj.mahal_cutoff_val)=0;

