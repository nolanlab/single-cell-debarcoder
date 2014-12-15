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
        cofactored_xt
        cofactored_xl
        mahal
        seprange
        clust_size

    end
    
    methods
        
        function obj = scd(key_filename)
            %constructor parses barcode key and sets default parameters
            
            %set defaults
            obj.default_cofactor=10;
            obj.sep_cutoff=0.3;
            obj.mahal_cutoff_val=30;
            
            axticks=load('axticks.mat');
            obj.xtl=axticks.xtl;
            obj.raw_xl=[-10, 10000];
            obj.raw_xt=axticks.xt;
            obj.default_xl=bmtrans(obj.raw_xl,obj.default_cofactor);
            obj.default_xt=bmtrans(obj.raw_xt,obj.default_cofactor);
            
            
            [pathstr, name, ext]=fileparts(key_filename);
            
            if ~strcmpi(ext,'.csv')
                error('Barcode key must be a csv file.')
            end
            
            %if can find file ...
            
            obj.key_filename=name;
            
            y=importdata(key_filename);
            obj.masses=cellstr(num2str(y.data(1,:)'));
            
            obj.masses=cellstr(num2str(y.data(1,:)'));
            obj.num_masses=length(obj.masses);
            obj.wellLabels=y.textdata(2:end);
            obj.num_codes=length(obj.wellLabels);
            
            obj.key=y.data(2:end,:);
            
            obj.well_yield=zeros(obj.num_codes,1);
            
            %add something here so that if fcs file already loaded,
            %it clears and updates
            
            % if don't recofactor
            obj.cofactored_xt=repmat(bmtrans(obj.raw_xt,obj.default_cofactor),[1 obj.num_masses]);
            obj.cofactored_xl=obj.cofactored_xt([1 end],:);
            obj.cofactors=repmat(obj.default_cofactor,[1 obj.num_masses]);
            
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
            
            %remove any data specific to previous fcs file
            obj.bcs=[];
            obj.cofactored_bcs=[];
            obj.well_yield=[];
            obj.bc_cols=[];
            obj.sample_ratio=[];
            obj.normbcs=[];
            obj.deltas=[];
            obj.bcind=[];
            obj.mahal=[];
            obj.clust_size=[];
            
        end
        
        function obj = find_bc_cols_by_mass(obj)
            
            obj.bc_cols=zeros(1,obj.num_masses);
            for i=1:obj.num_masses
                col_i=find(~cellfun(@isempty,regexp(obj.c,obj.masses(i))));
                if ~isempty(col_i) && length(col_i)==1
                    obj.bc_cols(i)=col_i;
                else
                    error('not all barcode channels found')
                end
            end
            
        end
        
        function obj = load_bcs(obj, sample_size)
        % extract and transform barcode columns from the fcs file based on the barcode key   
            if isempty(obj.x)
                error('An fcs file must be opened before loading BCs.')
            end
            
            num_cells=size(obj.x,1);
            
            if nargin==2 && num_cells>sample_size %sample
                obj.bcs=obj.x(randsample(num_cells,sample_size),obj.bc_cols);
                obj.sample_ratio=num_cells/sample_size;
            else %use all events
                obj.bcs=obj.x(:,obj.bc_cols);  
                obj.sample_ratio=1;
            end
            
            %transform
%             for i=1:obj.num_masses
%                 obj.bcs(:,i)=bmtrans(obj.bcs(:,i),obj.cofactors(i));
%             end
            obj.bcs=bmtrans(obj.bcs,obj.default_cofactor);
     
            %default, will change if recofactor
            obj.cofactored_bcs=obj.bcs;
        end
        
        function obj = normalize_bcs(obj,fieldname)
            % rescale the data in obj.fieldname for each parameter: note that this step
            % assumes that every parameter has both a positive and negative population
            % and when this assumption fails the rescaling can lead to wrong barcode
            % assignment 
            %fieldname should be 'bcs' or 'cofactored_bcs'
            
            eval(['data = obj.' fieldname ';'])
            
            percs=prctile(data,[1 99]);
            ranges=diff(percs);
            diffs=bsxfun(@minus,data,percs(1,:));
            obj.normbcs=bsxfun(@rdivide,diffs,ranges);
 
        end
        
        function obj = normalize_by_pop(obj,fieldname)
            %rescales the transformed bcs for each population using preliminary
            %assignments
            
            if nargin < 2
                fieldname = 'bcs';
            end
            
            eval(['data = obj.' fieldname ';'])
            
            bc_num_thresh=1;
            normed_bcs=zeros(size(data));
            for i=1:obj.num_codes
                inbc=obj.bcind==i;
                if nnz(inbc)>bc_num_thresh
                    pos_bcs=data(inbc,obj.key(i,:)==1);
                    %         norm_val=median(pos_bcs(:));
                    norm_val=prctile(pos_bcs(:),95);
                    normed_bcs(inbc,:)=data(inbc,:)/norm_val;
                    
                    %         for j=pos_inds
                    %            normed_bcs(inbc,j)=handles.bcs(inbc,j)/prctile(handles.bcs(inbc,j),95);
                    %         end
                    %         neg_inds=find(handles.key(i,:)==0);
                    %         for j=neg_inds
                    %            normed_bcs(inbc,j)=handles.bcs(inbc,j)-median(handles.bcs(inbc,j));
                    %         end
                end
                
            end
            
            obj.normbcs=normed_bcs;
            
        end
        
        function obj=compute_debarcoding(obj,fieldname)
            % for each event of normalized barcode intensities, assign that event to a
            % barcode and calculate the barcode separation
            
            %default to using normbcs
            if nargin<2
                fieldname='normbcs';
            end
            
            eval(['data = obj.' fieldname ';'])
            
            cutoff=0; %this was used to prevent large negative values from appearing to have sufficient separation from values near zero (not
            %needed without the +/-100 routine)
            N=size(obj.bcs,1);
            indlist=(1:N)';
            
            if length(unique(sum(obj.key,2)))==1
                % doublet-filtering code: look at top k barcodes, rather than barcodes
                % above largest separation
                
                [sorted,ix]=sort(data,2,'descend'); %barcode intensities ordered within each event
                
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
                
                [sorted,ix]=sort(data,2,'ascend'); %barcode intensities ordered within each event
                
                seps=diff(sorted,1,2); %barcode separations between every consecutive ordered barcode within each event
                
                [~,locs]=sort(seps,2,'descend');  %locs are columns in ix of bc level that is on lower side of largest gap, e.g., if locs is 5, the largest bc separation is between barcode ix(5) and ix(6)
                
                betws=ix(sub2ind(size(ix),indlist,locs(:,1)+1));  %columns of lowest barcode that is above the largest separation in each event
                lowests=data(sub2ind(size(obj.bcs),indlist,betws));  %normalized transformed values of lowest 'positive' BC
                
                betws=ix(sub2ind(size(ix),indlist,locs(:,1))); % columns of highest barcode that is below the largest separation in each event
                nextlowests=data(sub2ind(size(obj.bcs),indlist,betws)); %normalized transformed values of highest 'negative' BC
                
                toolow=find(obj.bcs(indsabove)<cutoff);  %these aren't high enough to count.  go to next-biggest-sep. using actual bcs here, not normalized
                
                betws=ix(sub2ind(size(ix),toolow,locs(toolow,2)+1));
                inds=sub2ind(size(obj.bcs),toolow,betws);
                lowests_next=data(inds);
                highernow=obj.bcs(inds)>cutoff;  %again using actualy bcs, not normalized, to check against cutoff
                %might still need to account for when the largest sep is high ...  can
                %first try eliminating these just with illegal barcodes
                lowests(toolow(highernow))=lowests_next(highernow);
                lowests(toolow(~highernow))=nan;  %if the second try didn't find one above the cutoff, set to nan
                
                %adding in the replaced bcs
                betws=ix(sub2ind(size(ix),toolow,locs(toolow,2)));
                inds=sub2ind(size(obj.bcs),toolow,betws);
                
                modifiednextlowests2=data(inds);
                nextlowests(toolow(highernow))=nextlowests2(highernow);
                nextlowests(toolow(~highernow))=nan;
                
                deltas=lowests-nextlowests;  %separation between 'positive' and 'negative' barcodes for each cell
                
            end
            
            obj.deltas=deltas;
            
            % assign binary barcodes to each cell
            code_assign=false(N,obj.num_masses);
            for j=1:obj.num_masses
                code_assign(:,j)=data(:,j) >= lowests;
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
        
        function obj=calculate_cofactors(obj)
            % determine a cofactor for each bc channel by pooling the negative barcodes
            % for that channel across the populations
            
            temp_bcind=obj.bcind;
%             temp_bcind(obj.deltas<obj.sep_cutoff)=0;
            N=length(temp_bcind);
            
            neg_bcs=cell(1,obj.num_masses);
            bc_list=1:obj.num_codes;
            neg_cofactor=zeros(1,obj.num_masses);
            for i=1:obj.num_masses
                neg_list=bc_list(obj.key(:,i)==0);
                neg_cells=false(N,1);
                for j=neg_list
                    neg_cells=neg_cells | temp_bcind==j;
                end

                neg_bcs{i}=obj.bcs(neg_cells,i); %this was already transformed using default cofactor
                neg_cofactor(i)=obj.default_cofactor*sinh(prctile(neg_bcs{i},50)); %untransformed to raw data val
            end
            
            if any(isnan(neg_cofactor))
                warndlg('Check your barcode key. You may have included a barcode metal that is constant across all occupied samples.')
                neg_cofactor(isnan(neg_cofactor))=5;
            end
            
            neg_cofactor(neg_cofactor<5)=5; %5 is default minimum
            neg_cofactor(neg_cofactor>100)=100; %100 is default maximum ... maybe should lower
            obj.cofactors=neg_cofactor;
        end      
        
        function obj=recofactor(obj)
            %retransform cofactored_bcs and norm_vals from default cofactoring to variable
            %cofactoring
            
%             cofactored_bcs=zeros(size(obj.bcs)); %already exists from
%             when loaded
            
            for i=1:obj.num_masses
                cofactored_bcs(:,i)=asinh(obj.default_cofactor*sinh(obj.bcs(:,i))/obj.cofactors(i));
                obj.cofactored_xt(:,i)=bmtrans(obj.raw_xt,obj.cofactors(i));
            end
            obj.cofactored_xl=obj.cofactored_xt([1 end],:);
            
            obj.cofactored_bcs=cofactored_bcs;

        end
        
        function obj = compute_mahal(obj)
            % computes the mahalanobis distances of all events given the current
            % separation cutoff, right now uses bcs but maybe should use
            % cofactored or norm?
            
            obj.mahal=zeros(size(obj.deltas));
            for i=1:obj.num_codes
                in_bc=(obj.bcind==i) & (obj.deltas > obj.sep_cutoff);
                bci=obj.bcs(in_bc,:);
                if size(bci,1)>obj.num_codes
                    obj.mahal(in_bc)=mahal(bci,bci);
                end
                obj.well_yield(i)=obj.sample_ratio*nnz(in_bc & obj.mahal<obj.mahal_cutoff_val);
            end

        end
        
        function obj=compute_well_abundances(obj)
            % compute well abundances
            
            numseps=20;
            minsep=0;
            maxsep=1;
            obj.seprange=linspace(minsep,maxsep,numseps);
            
            obj.clust_size=zeros(numseps,obj.num_codes);
            
            for i=1:numseps
                for j=1:length(obj.wellLabels)
                    obj.clust_size(i,j) = nnz(obj.bcind==j & (obj.deltas > obj.seprange(i)));
                end
            end
            
        end
        
        function write_bc_fcs_files(obj,outdir,basename)
            %write an fcs for each barcode population, and a file of
            %unassigned events. save in the directory
            %outdir.
            
            
            % write an fcs file for each barcode
            not_inawell=true(size(obj.bcind));
            for i=1:obj.num_codes
                thiswell_bin = (obj.bcind==i) & (obj.mahal<obj.mahal_cutoff_val) & (obj.deltas > obj.sep_cutoff);
                not_inawell(thiswell_bin)=false; %cells in this well removed from unassigned_binary
                data=obj.x(thiswell_bin,:);
                if ~isempty(data)
                    fca_writefcs([outdir filesep basename '_' obj.wellLabels{i} '.fcs'],data,obj.m,obj.c)
                end
            end
            
            % write fcs file of unassigned events
            data=obj.x(not_inawell,:);
            fca_writefcs([outdir filesep basename '_unassigned.fcs'],data,obj.m,obj.c)            
           
        end
        
    end
    
    methods(Static)
        
        function y=bmtrans(x,c)
            %asinh transform with cofactor c
            num_cols=size(x,2);
            if length(c)==1
                c=repmat(c,[1 num_cols]);
            end
            y=zeros(size(x));
            for i=1:num_cols
                y(:,i)=asinh(1/c(i)*x(:,i));
            end
        end
        
        
            
    end
    
    
end