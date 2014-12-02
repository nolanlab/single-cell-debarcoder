classdef scd
    
    properties
        sep_cutoff
        mahal_cutoff_val
        default_cofactor
        xtl
        raw_xl
        raw_xt
        default_xl
        default_xt
        key_filename
        masses
        num_masses
        wellLabels
        num_codes
        key
        well_yield
        x
        current_files
        c
        m
        bc_cols
        bcs
        sample_ratio
        normbcs
        deltas
        bcind
        cofactored_bcs
        cofactors
        cotactored_xt
        cofactored_xl
        mahal
        seprange
        clust_size
        ax
        lines
        leg
    end
    
    methods
        
        function obj = scd(key_filename)
            
            [pathstr, name, ext]=fileparts(key_filename);
            
            if ~strcmpi(ext,'.csv')
                error('Barcode key must be a csv file.')
            end
            
            %if can find file ...
            
            obj.key_filename=name;
            
            x=importdata(key_filename);
            obj.masses=cellstr(num2str(x.data(1,:)'));
            
            obj.masses=cellstr(num2str(x.data(1,:)'));
            obj.num_masses=length(obj.masses);
            obj.wellLabels=x.textdata(2:end);
            obj.num_codes=length(obj.wellLabels);
            
            obj.key=x.data(2:end,:);
            
            obj.well_yield=zeros(obj.num_codes,1);
            
            %add something here so that if fcs file already loaded,
            %it clears and updates
            
            
        end
        
        function obj=load_fcs_files(obj,filenames)
            
            if iscell(filenames) %more than 1 file selected, which means concatenation
                num_files=length(filenames);
                z=cell(1,num_files);
                for i=1:num_files
                    [z{i},h]=fca_readfcs(filenames{i});
                end
                obj.x=cat(1,z{:});
                
                
            else
                [obj.x,h]=fca_readfcs(filenames);
                
            end
            
            obj.c={h.par.name};
            obj.m={h.par.name2};
            
        end
        
        function obj=compute_debarcoding(obj)
            % for each event of normalized barcode intensities, assign that event to a
            % barcode and calculate the barcode separation
            
            cutoff=0; %this was used to prevent large negative values from appearing to have sufficient separation from values near zero (not
            %needed without the +/-100 routine)
            N=size(obj.bcs,1);
            indlist=(1:N)';
            
            if length(unique(sum(obj.key,2)))==1
                % doublet-filtering code: look at top k barcodes, rather than barcodes
                % above largest separation
                
                [sorted,ix]=sort(obj.normbcs,2,'descend'); %barcode intensities ordered within each event
                
                numdf=sum(obj.key(1,:));  %number of expected positive barcode intensities
                
                lowests=sorted(:,numdf); %the value of the lowest 'positive' BC for each cell
                
                %get rid of cells whose 'positive' barcodes are still very low (not
                %needed without the +/-100 routine)
                inds=sub2ind(size(ix),indlist,ix(:,numdf));
                toolow=obj.bcs(inds)<cutoff; %using bcs, not normbcs
                lowests(toolow)=nan;
                
                deltas=sorted(:,numdf)-sorted(:,numdf+1);  %separation between 'positive' and 'negative' barcodes for each cell
                
            else
                % non-constant number of '1's in code, so find largest separation within each event to assign 'positive' and 'negative' channels
                
                [sorted,ix]=sort(obj.normbcs,2,'ascend'); %barcode intensities ordered within each event
                
                seps=diff(sorted,1,2); %barcode separations between every consecutive ordered barcode within each event
                
                [~,locs]=sort(seps,2,'descend');  %locs are columns in ix of bc level that is on lower side of largest gap, e.g., if locs is 5, the largest bc separation is between barcode ix(5) and ix(6)
                
                betws=ix(sub2ind(size(ix),indlist,locs(:,1)+1));  %columns of lowest barcode that is above the largest separation in each event
                lowests=obj.normbcs(sub2ind(size(obj.bcs),indlist,betws));  %normalized transformed values of lowest 'positive' BC
                
                betws=ix(sub2ind(size(ix),indlist,locs(:,1))); % columns of highest barcode that is below the largest separation in each event
                nextlowests=obj.normbcs(sub2ind(size(obj.bcs),indlist,betws)); %normalized transformed values of highest 'negative' BC
                
                toolow=find(obj.bcs(indsabove)<cutoff);  %these aren't high enough to count.  go to next-biggest-sep. using actual bcs here, not normalized
                
                betws=ix(sub2ind(size(ix),toolow,locs(toolow,2)+1));
                inds=sub2ind(size(obj.bcs),toolow,betws);
                lowests_next=obj.normbcs(inds);
                highernow=obj.bcs(inds)>cutoff;  %again using actualy bcs, not normalized, to check against cutoff
                %might still need to account for when the largest sep is high ...  can
                %first try eliminating these just with illegal barcodes
                lowests(toolow(highernow))=lowests_next(highernow);
                lowests(toolow(~highernow))=nan;  %if the second try didn't find one above the cutoff, set to nan
                
                %adding in the replaced bcs
                betws=ix(sub2ind(size(ix),toolow,locs(toolow,2)));
                inds=sub2ind(size(obj.bcs),toolow,betws);
                
                modifiednextlowests2=obj.normbcs(inds);
                nextlowests(toolow(highernow))=nextlowests2(highernow);
                nextlowests(toolow(~highernow))=nan;
                
                deltas=lowests-nextlowests;  %separation between 'positive' and 'negative' barcodes for each cell
                
            end
            
            obj.deltas=deltas;
            
            % assign binary barcodes to each cell
            code_assign=false(N,obj.num_masses);
            for j=1:obj.num_masses
                code_assign(:,j)=obj.normbcs(:,j) >= lowests;
            end
            
            % assign barcode ID (1:num_barcodes) to each cell
            obj.bcind=zeros(N,1);
            num_cells=size(obj.bcs,1);
            for i=1:obj.num_codes
                clust_inds=true(num_cells,1);
                for j=1:obj.num_masses
                    clust_inds=clust_inds & (code_assign(:,j)==obj.key(i,j));
                end
                obj.bcind(clust_inds)=i;
            end
        end
    end
    
end