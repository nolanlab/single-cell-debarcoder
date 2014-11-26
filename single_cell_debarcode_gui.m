function varargout = single_cell_debarcode_gui(varargin)
% SINGLE_CELL_DEBARCODE_GUI MATLAB code for single_cell_debarcode_gui.fig
%      SINGLE_CELL_DEBARCODE_GUI, by itself, creates a new SINGLE_CELL_DEBARCODE_GUI or raises the existing
%      singleton*.
%
%      H = SINGLE_CELL_DEBARCODE_GUI returns the handle to a new SINGLE_CELL_DEBARCODE_GUI or the handle to
%      the existing singleton*.
%
%      SINGLE_CELL_DEBARCODE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SINGLE_CELL_DEBARCODE_GUI.M with the given input arguments.
%
%      SINGLE_CELL_DEBARCODE_GUI('Property','Value',...) creates a new SINGLE_CELL_DEBARCODE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before single_cell_debarcode_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to single_cell_debarcode_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help single_cell_debarcode_gui

% Last Modified by GUIDE v2.5 22-May-2014 14:48:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @single_cell_debarcode_gui_OpeningFcn, ...
    'gui_OutputFcn',  @single_cell_debarcode_gui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before single_cell_debarcode_gui is made visible.
function single_cell_debarcode_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to single_cell_debarcode_gui (see VARARGIN)

% Choose default command line output for single_cell_debarcode_gui
handles.output = hObject;

% color order
co=[     0    0.5000    0.4000
    0.7500         0    0.7500
    1.0000    0.8000         0
         0         0    0.7500
    1.0000    0.3750         0
         0    0.7500    1.0000
    0.4000    0.4000    0.4000];

% parent figure properties
handles.parent=gcf;
set(handles.parent,'WindowStyle','normal','Name','Single Cell Debarcoder',...
    'DefaultAxesColorOrder',co)

% setup default initial view
set(handles.well_popup,'visible','off')
set(handles.pop_fwd,'visible','off')
set(handles.pop_back,'visible','off')
set(handles.biax_panel,'visible','off')
set(handles.color_panel,'visible','off')
set(handles.plottype,'SelectionChangeFcn',{@plot_changefcn,handles})
set(handles.plottype,'SelectedObject',handles.separation_plot)
set(handles.color_panel,'SelectionChangeFcn',{@color_changefcn,handles})
set(handles.color_panel,'SelectedObject',handles.color_mahal)
set(handles.plot_button,'visible','off')

% initialize parameters
set(handles.delta_text,'string','0.3')
handles.sep_cutoff=0.3;

set(handles.mahal_cutoff,'string','30')
handles.mahal_cutoff_val=30;

handles.default_cofactor=5;

% axis limits and ticks
axticks=load('axticks.mat');
handles.xtl=axticks.xtl;
handles.raw_xl=[-10, 10000];
handles.raw_xt=axticks.xt;
handles.default_xl=bmtrans(handles.raw_xl,handles.default_cofactor);
handles.default_xt=bmtrans(axticks.xt,handles.default_cofactor);

% remove unwanted toolbar options
set(hObject,'toolbar','figure');
hu=findall(hObject,'type','uitoolbar');
ch=findall(hu);
chbin=false(length(ch),1);
tags={'FigureToolBar','Exploration.DataCursor','Exploration.Pan','Exploration.ZoomOut','Exploration.ZoomIn'};
for i=1:length(ch)
    if ~any(strcmp(get(ch(i),'tag'),tags))
        chbin(i)=true;
    end
end
delete(ch(chbin))

% update handles
guidata(hObject, handles)

% load barcode key
bc_button_Callback(hObject, eventdata, handles);


% --- Outputs from this function are returned to the command line.
function varargout = single_cell_debarcode_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function color_changefcn(hObject,eventdata,handles)
% executes upon change in selection of the 'Color Plot By' radio button.
% This sets the coloring parameter to either Mahalanobis distance or
% Barcode Separation Threshold and re-plots the data.
handles=guidata(handles.parent);
plot_button_Callback(hObject,eventdata,handles)

function plot_changefcn(hObject,eventdata,handles)
% executes upon change in selection of the 'Plot Type' radio button.
% Switches the view to the selected plot type and re-plots the data.
selectedobj = get(handles.plottype,'SelectedObject');
handles=guidata(handles.parent);
switch get(selectedobj,'string')
    case 'Separation'
        set(handles.well_popup,'visible','off')
        set(handles.pop_fwd,'visible','off')
        set(handles.pop_back,'visible','off')
        set(handles.biax_panel,'visible','off')
        set(handles.color_panel,'visible','off')        
    case 'Event'
        set(handles.well_popup,'visible','on')
        set(handles.pop_fwd,'visible','on')
        set(handles.pop_back,'visible','on')
        set(handles.biax_panel,'visible','off')
        set(handles.color_panel,'visible','off')
        plot_well_yields(handles)
    case 'Single Biaxial'
        set(handles.well_popup,'visible','on')
        set(handles.pop_fwd,'visible','on')
        set(handles.pop_back,'visible','on')
        set(handles.biax_panel,'visible','on')
        set(handles.color_panel,'visible','on')
        plot_well_yields(handles)
    case 'All Barcode Biaxials'
        set(handles.well_popup,'visible','on')
        set(handles.pop_fwd,'visible','on')
        set(handles.pop_back,'visible','on')
        set(handles.biax_panel,'visible','off')
        set(handles.color_panel,'visible','on')
        plot_well_yields(handles)
end
plot_button_Callback(hObject, eventdata, handles)


function plot_button_Callback(hObject, eventdata, handles)
% creates the plot of the type selected by the 'Plot Type' radio button and
% the parameters selected in the 'Plot Options' panel. The plot button has
% been removed and instead this is set to execute automatically upon any change to plot
% type or plot options.

selectedobj = get(handles.plottype,'SelectedObject');

wellnum=get(handles.well_popup,'value');
thiswell_bin = (handles.bcind==wellnum) & (handles.mahal<handles.mahal_cutoff_val) & (handles.deltas > handles.sep_cutoff);

if ~isempty(wellnum)
    
    set(handles.parent,'pointer','watch')
    drawnow;
    switch get(selectedobj,'string')
        case 'Separation'
            separation_plot(handles);
        case 'Event'
            
            oldax=get(handles.ax_panel,'children');
            delete(oldax)
            
            %raw bcs
            ax=zeros(1,2);
            ax(1)=subplot(2,1,1,'parent',handles.ax_panel);
            hold on
            
            thiswell=handles.bcs(thiswell_bin,:);
            if ~isempty(thiswell)
                s=1:size(thiswell,1);
                plot3(s,thiswell,rand(size(s))-1,'.')
                hold on              
                set(ax(1),'Box','off','xtick',[],'Layer','top')
                set(ax(1),'ytick',handles.default_xt,'yticklabel',handles.xtl)
                set(get(ax(1),'Ylabel'),'String','untransformed values','fontsize',12)
                set(get(ax(1),'Ylabel'),'String','Barcode intensities','fontsize',12)
                xlabel('Event','fontsize',12)
                legend(handles.leg,'location','northeastoutside')    
            end
            title(num2str(handles.key(wellnum,:)),'fontweight','bold','fontsize',12);
            
            %norm bcs
            ax(2)=subplot(2,1,2,'parent',handles.ax_panel);
            hold on
            
            thiswell=handles.normbcs(thiswell_bin,:);           
            if ~isempty(thiswell)
                s=1:size(thiswell,1);
                plot3(s,thiswell,rand(size(s)),'.');
                set(ax(2),'xtick',[])
                set(get(ax(2),'Ylabel'),'String','Rescaled values','fontsize',12)
                xlabel('Event','fontsize',12)
                legend(handles.leg,'location','northeastoutside')
            end
            
            hlink=linkprop(ax,'xLim');
            key = 'graphics_linkprop';
            setappdata(ax(1),key,hlink);
                        
        case 'Single Biaxial'
            oldax=get(handles.ax_panel,'children');
            delete(oldax)
            handles.ax=axes('parent',handles.ax_panel);
            xcol=get(handles.x_popup,'value');
            ycol=get(handles.y_popup,'value');
            cm=jet(64);
            hold on
            
            thiswell=handles.cofactored_bcs(thiswell_bin,[xcol ycol]);
            
            mdists=handles.mahal(thiswell_bin);
            
            seps=handles.deltas(thiswell_bin);
            if ~isempty(thiswell)
                
                if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                    scatter(thiswell(:,1),thiswell(:,2),4,mdists)
                    set(gca,'clim',[0 handles.mahal_cutoff_val])
                    colormap(flipud(cm))
                else
                    scatter(thiswell(:,1),thiswell(:,2),4,seps)
                    set(gca,'clim',[handles.sep_cutoff 1])
                    
                    colormap(cm)
                end
                
                cb=colorbar;
                if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                    set(get(cb,'ylabel'),'string','Mahalanobis Distance','fontsize',14)
                else
                    set(get(cb,'ylabel'),'string','Separation','fontsize',14)
                end
                
            end
            
            set(gca,'xlim',handles.cofactored_xl(:,xcol),'ylim',handles.cofactored_xl(:,ycol),...
                'xtick',handles.cofactored_xt(:,xcol),'xticklabel',handles.xtl,...
                'ytick',handles.cofactored_xt(:,ycol),'yticklabel',handles.xtl)
            
            title(num2str(handles.key(wellnum,:)),'fontweight','bold','fontsize',12);
            
        case 'All Barcode Biaxials'
            oldax=get(handles.ax_panel,'children');
            delete(oldax)
            n=size(handles.bcs,2);
            
            cm=jet(64);
            
            Ind=1;
            for i=0:n
                for j=0:n
                    
                    if j==0 && i~=n
                        ax=subplot(n+1,n+1,Ind,'Parent',handles.ax_panel);
                        set(ax,'visible','off')
                        text(0.5,0.5,handles.masses{i+1},'horizontalalignment','center','verticalalignment','middle')
                    elseif j<=i && i<n
                        ax=subplot(n+1,n+1,Ind,'Parent',handles.ax_panel);
                        hold on

                        thiswell=handles.cofactored_bcs(thiswell_bin,[j i+1]);

                        mdists=handles.mahal(thiswell_bin);
                        seps=handles.deltas(thiswell_bin);
                        
                        if ~isempty(thiswell)
                            if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                                scatter(thiswell(:,1),thiswell(:,2),2,mdists)
                                set(gca,'clim',[0 handles.mahal_cutoff_val])
                            else
                                scatter(thiswell(:,1),thiswell(:,2),2,seps)
                                set(gca,'clim',[handles.sep_cutoff 1])
                            end
                        end
                        
                        set(ax,'xlim',handles.cofactored_xl(:,j),'ylim',handles.cofactored_xl(:,i+1),'xtick',[],'ytick',[])
                    elseif j>0 && i==j-1 
                        ax=subplot(n+1,n+1,Ind,'Parent',handles.ax_panel);
                        if nnz(handles.bcind==wellnum) ~= 0
                            [binsize,binloc]=hist(handles.cofactored_bcs(thiswell_bin,j),100);
                            bar(binloc,binsize,'edgecolor','none','facecolor',[0 0.5 0])
                        end
                        set(ax,'xlim',handles.cofactored_xl(:,j),'xtick',[],'ytick',[])
                    elseif i==n && j>0
                        ax=subplot(n+1,n+1,Ind,'Parent',handles.ax_panel);
                        set(ax,'visible','off')
                        text(0.5,0.5,handles.masses{j},'horizontalalignment','center','verticalalignment','middle')
                    end
                    Ind=Ind+1;
                end
            end
            
            bigax=axes('parent',handles.ax_panel);
            if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                colormap(flipud(cm))
                set(bigax,'clim',[0 handles.mahal_cutoff_val])
                cb=colorbar('position',[0.85 0.5 0.027 0.4]);
                set(get(cb,'ylabel'),'string','Mahalanobis Distance','fontsize',14)
            else
                colormap(cm)
                set(bigax,'clim',[handles.sep_cutoff 1])
                cb=colorbar('position',[0.85 0.5 0.027 0.4]);
                set(get(cb,'ylabel'),'string','Separation','fontsize',14)
            end
            set(bigax,'visible','off')
    end
    
    drawnow
    set(handles.parent,'pointer','arrow')
    drawnow
    
else
    msgbox('No barcode with that label exists!')
end

function well_text_Callback(hObject, eventdata, handles)
% hObject    handle to well_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of well_text as text
%        str2double(get(hObject,'String')) returns contents of well_text as a double


% --- Executes during object creation, after setting all properties.
function well_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to well_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles=compute_debarcoding(handles)
% for each event of normalized barcode intensities, assign that event to a
% barcode and calculate the barcode separation

cutoff=0; %this was used to prevent large negative values from appearing to have sufficient separation from values near zero (not
    %needed without the +/-100 routine)
N=size(handles.bcs,1);
indlist=(1:N)';

if length(unique(sum(handles.key,2)))==1 
% doublet-filtering code: look at top k barcodes, rather than barcodes
% above largest separation
    
    [sorted,ix]=sort(handles.normbcs,2,'descend'); %barcode intensities ordered within each event
    
    numdf=sum(handles.key(1,:));  %number of expected positive barcode intensities 
    
    lowests=sorted(:,numdf); %the value of the lowest 'positive' BC for each cell
    
    %get rid of cells whose 'positive' barcodes are still very low (not
    %needed without the +/-100 routine)
    inds=sub2ind(size(ix),indlist,ix(:,numdf));
    toolow=handles.bcs(inds)<cutoff; %using bcs, not normbcs
    lowests(toolow)=nan;    
    
    deltas=sorted(:,numdf)-sorted(:,numdf+1);  %separation between 'positive' and 'negative' barcodes for each cell
    
else
% non-constant number of '1's in code, so find largest separation within each event to assign 'positive' and 'negative' channels
    
    [sorted,ix]=sort(handles.normbcs,2,'ascend'); %barcode intensities ordered within each event
    
    seps=diff(sorted,1,2); %barcode separations between every consecutive ordered barcode within each event
    
    [~,locs]=sort(seps,2,'descend');  %locs are columns in ix of bc level that is on lower side of largest gap, e.g., if locs is 5, the largest bc separation is between barcode ix(5) and ix(6)
    
    betws=ix(sub2ind(size(ix),indlist,locs(:,1)+1));  %columns of lowest barcode that is above the largest separation in each event
    lowests=handles.normbcs(sub2ind(size(handles.bcs),indlist,betws));  %normalized transformed values of lowest 'positive' BC
    
    betws=ix(sub2ind(size(ix),indlist,locs(:,1))); % columns of highest barcode that is below the largest separation in each event
    nextlowests=handles.normbcs(sub2ind(size(handles.bcs),indlist,betws)); %normalized transformed values of highest 'negative' BC
    
    toolow=find(handles.bcs(indsabove)<cutoff);  %these aren't high enough to count.  go to next-biggest-sep. using actual bcs here, not normalized
    
    betws=ix(sub2ind(size(ix),toolow,locs(toolow,2)+1));
    inds=sub2ind(size(handles.bcs),toolow,betws);
    lowests_next=handles.normbcs(inds);
    highernow=handles.bcs(inds)>cutoff;  %again using actualy bcs, not normalized, to check against cutoff
    %might still need to account for when the largest sep is high ...  can
    %first try eliminating these just with illegal barcodes
    lowests(toolow(highernow))=lowests_next(highernow);
    lowests(toolow(~highernow))=nan;  %if the second try didn't find one above the cutoff, set to nan
    
    %adding in the replaced bcs
    betws=ix(sub2ind(size(ix),toolow,locs(toolow,2)));
    inds=sub2ind(size(handles.bcs),toolow,betws);

    modifiednextlowests2=handles.normbcs(inds);
    nextlowests(toolow(highernow))=nextlowests2(highernow);
    nextlowests(toolow(~highernow))=nan;
    
    deltas=lowests-nextlowests;  %separation between 'positive' and 'negative' barcodes for each cell
    
end

handles.deltas=deltas;

% assign binary barcodes to each cell 
code_assign=false(N,handles.num_masses);  
for j=1:handles.num_masses
    code_assign(:,j)=handles.normbcs(:,j) >= lowests;
end

% assign barcode ID (1:num_barcodes) to each cell
handles.bcind=zeros(N,1);
num_cells=size(handles.bcs,1);
for i=1:handles.num_codes
    clust_inds=true(num_cells,1);
    for j=1:handles.num_masses
        clust_inds=clust_inds & (code_assign(:,j)==handles.key(i,j));
    end
    handles.bcind(clust_inds)=i;
end

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


FileName=get(handles.save_text,'string');

sep_cutoff=handles.sep_cutoff;

PathName = uigetdir;

if PathName ~= 0
    
    set(handles.parent,'pointer','watch')
    drawnow
    
    % recompute handles.bcind using all (not just sampled) cells    
    handles.bcs=zeros(size(handles.x,1),handles.num_masses);
    for i=1:handles.num_masses
    handles.bcs(:,i)=bmtrans(handles.x(:,handles.bc_cols(i)),handles.cofactors(i)); 
    end
    
    handles.normbcs=normalize_bcs(handles.bcs);

    handles=compute_debarcoding(handles);
    
    handles=compute_mahal(handles);
    
    % write an fcs file for each barcode
    not_inawell=true(size(handles.bcind));
    for i=1:length(handles.wellLabels)    
        thiswell_bin = (handles.bcind==i) & (handles.mahal<handles.mahal_cutoff_val) & (handles.deltas > sep_cutoff);        
        not_inawell(thiswell_bin)=false; %cells in this well removed from unassigned_binary       
        data=handles.x(thiswell_bin,:);
        if ~isempty(data)
            fca_writefcs([PathName filesep FileName '_' handles.wellLabels{i} '.fcs'],data,handles.m,handles.c)
        end
    end

    % write fcs file of unassigned events
    data=handles.x(not_inawell,:);
    fca_writefcs([PathName filesep FileName '_unassigned.fcs'],data,handles.m,handles.c)

    % write log file of debarcoding process
    fid=fopen([PathName filesep 'Debarcode_Parameters.txt'],'w');
    fprintf(fid,'%s\n',datestr(now));
    fprintf(fid,'%s\n','Input files:');
    for i=1:length(handles.current_files)
        fprintf(fid,'\t%s\n',handles.current_files{i});
    end
    fprintf(fid,'%s\n',['Barcode Key: ' handles.key_filename]);
    fprintf(fid,'%s\n',['Mahalanobis cutoff: ' get(handles.mahal_cutoff,'string')]);
    fprintf(fid,'%s\n',['Separation cutoff: ' get(handles.delta_text,'string')]);
    fprintf(fid,'%s\n',['Output saved to: ' PathName filesep]);
    fclose(fid);
    
    set(handles.parent,'pointer','arrow')
    drawnow
else
    return
end

function normed_bcs = normalize_bcs(raw_bcs)
% rescale the transformed bcs for each parameter: note that this step
% assumes that every parameter has both a positive and negative population
% and when this assumption fails the rescaling can lead to wrong barcode
% assignment

percs=prctile(raw_bcs,[1 99]);
ranges=diff(percs);
deltas=bsxfun(@minus,raw_bcs,percs(1,:));
normed_bcs=bsxfun(@rdivide,deltas,ranges);



% --- Executes on button press in folder_button.
function folder_button_Callback(hObject, eventdata, handles)
% read in an fcs file. if multiple files are selected, they are
% concatenated

[files,handles.pathname] = uigetfile('*.fcs','Select fcs files to debarcode','multiselect','on');

if handles.pathname ~= 0    
    set(handles.parent,'pointer','watch')
    drawnow
    
    if iscell(files) %more than 1 file selected, which means concatenation
        num_files=length(files);
        z=cell(1,num_files);
        for i=1:num_files
            [z{i},h]=fca_readfcs([handles.pathname filesep files{i}]);
        end
        handles.x=cat(1,z{:});
       
        formatstr=repmat('\n%s',[1 num_files]);
        str=sprintf(['Using files:' formatstr],files{:});        
        handles.current_files=files;
    else
        [handles.x,h]=fca_readfcs([handles.pathname filesep files]);
        
        str=sprintf('Using file:\n%s',files);        
        handles.current_files={files};
    end    
    
    handles.c={h.par.name};
    handles.m={h.par.name2};
    
    set(handles.file_text,'string',str);
    
    % find bc columns and do preliminary transformation
    handles=load_bc_data(hObject, eventdata,handles);
    
    guidata(hObject,handles)
    
else
    return
end

function handles=compute_well_abundances(handles)
% compute well abundances

numseps=20;
minsep=0;
maxsep=1;
handles.seprange=linspace(minsep,maxsep,numseps);

handles.clust_size=zeros(numseps,handles.num_codes);

for i=1:numseps
    for j=1:length(handles.wellLabels)
        handles.clust_size(i,j) = nnz(handles.bcind==j & (handles.deltas > handles.seprange(i)));
    end
end

function handles=separation_plot(handles)
% makes a histogram of total yield binned by barcode separation, and a plot of each well's yield as
% a function of separation cutoff

ch=get(handles.yield_panel,'children');
delete(ch)
yield_axis=axes('parent',handles.yield_panel);
[hi,xi]=hist(handles.deltas,100);
bar(yield_axis,xi,handles.sample_ratio*hi,'facecolor',[0 0.5 0.4])
set(get(yield_axis,'XLabel'),'String','Barcode separation','fontsize',12)
set(get(yield_axis,'YLabel'),'String','Event count','fontsize',12)
set(yield_axis,'xlim',[0 1],'box','off')

ch=get(handles.ax_panel,'children');
delete(ch)
handles.ax=axes('parent',handles.ax_panel);
set(handles.ax,'colororder',flipud(jet(20)))
hold on
handles.lines=plot(handles.ax,handles.seprange,handles.sample_ratio*handles.clust_size);
set(get(handles.ax,'XLabel'),'String','Barcode separation threshold','fontsize',12)
set(get(handles.ax,'YLabel'),'String','Event yield after debarcoding','fontsize',12)
set(handles.lines,'ButtonDownFcn',{@select_line,handles});


function select_line(button,eventdata,handles)
% displays on the separation plot of barcode label of any line selected by
% the mouse

ch=get(handles.ax,'children');
T=ch(strcmp('text',get(ch,'type')));
if ~isempty(T)
    delete(T)
end

ind=handles.lines==gco;
cp=get(handles.ax,'currentpoint');
text(cp(1,1),cp(1,2),handles.wellLabels{ind},'fontsize',12,'verticalalignment','baseline','edgecolor','k','backgroundcolor','w','interpreter','none')

function handles=compute_mahal(handles)
% computes the mahalanobis distances of all events given the current
% separation cutoff

handles.mahal=zeros(size(handles.deltas));
for i=1:handles.num_codes
    in_bc=(handles.bcind==i) & (handles.deltas > handles.sep_cutoff);
    bci=handles.bcs(in_bc,:);
    if size(bci,1)>handles.num_codes
        handles.mahal(in_bc)=mahal(bci,bci);
    end
    handles.well_yield(i)=handles.sample_ratio*nnz(in_bc & handles.mahal<handles.mahal_cutoff_val);
end

function handles=delta_text_Callback(hObject, eventdata, handles)
% executes upon a change to the separation cutoff value

% update separation cutoff
handles.sep_cutoff=str2double(get(handles.delta_text,'string'));

% compute mahalanobis distances
handles=compute_mahal(handles);

% re-plots the data if a separation-cutoff dependent plot is being
% displayed
selectedobj=get(handles.plottype,'SelectedObject');
if ~strcmp(get(selectedobj,'string'),'Separation')
    plot_well_yields(handles)
    plot_button_Callback(hObject, eventdata, handles)
end

guidata(handles.parent,handles)


% --- Executes during object creation, after setting all properties.
function delta_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delta_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cutoff_text_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff_text as text
%        str2double(get(hObject,'String')) returns contents of cutoff_text as a double


% --- Executes during object creation, after setting all properties.
function cutoff_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in bc_button.
function bc_button_Callback(hObject, eventdata, handles) %"Change" button to select barcode key
% load a barcode key

[handles.key_filename, handles.bcpathname]=uigetfile({'*.csv','*.CSV'},'Select barcode key');

if handles.bcpathname ~= 0
    
    x=importdata([handles.bcpathname handles.key_filename]);
    handles.masses=cellstr(num2str(x.data(1,:)'));
    str=sprintf('Using key:\n%s',handles.key_filename);
    set(handles.key_text,'string',str);
    
    handles.masses=cellstr(num2str(x.data(1,:)'));
    handles.num_masses=length(handles.masses);
    handles.wellLabels=x.textdata(2:end);
    handles.num_codes=length(handles.wellLabels);
    
    set(handles.well_popup,'string',handles.wellLabels)
    
    handles.key=x.data(2:end,:);
    
    handles.well_yield=zeros(handles.num_codes,1);
       
    %if fcs file already loaded, update which are the bc cols and the bc data i
    if isfield(handles,'c')
        handles=load_bc_data(hObject, eventdata,handles);
    end
    guidata(hObject,handles)

else
    return
end

function handles=load_bc_data(hObject, eventdata, handles)
% extract barcode columns from the fcs file based on the barcode key

set(handles.parent,'pointer','watch')
drawnow

try
    handles.bc_cols=find_bc_cols_by_mass(handles.c,handles.masses);
catch err
    if strcmp('not all barcode channels found',err.message)
        set(handles.parent,'pointer','arrow')
        msgbox('Check your barcode key.')
        guidata(handles.parent,handles)
        return
    else
        set(handles.parent,'pointer','arrow')
        rethrow(err)
    end
end

%sample 100000 cells and use until save
num_cells=size(handles.x,1);
sample_size=100000;

if num_cells>sample_size
    handles.bcs=bmtrans(handles.x(randsample(num_cells,sample_size),handles.bc_cols),handles.default_cofactor);
    handles.sample_ratio=num_cells/sample_size;
else
    handles.bcs=bmtrans(handles.x(:,handles.bc_cols),handles.default_cofactor);  %matrix of each cell's bc channels, transformed
    handles.sample_ratio=1;
end


handles.normbcs=normalize_bcs(handles.bcs);

handles=compute_debarcoding(handles);

%% 20140904 -- main cofactor update 
handles=calculate_cofactors(handles);

handles=recofactor(handles);

handles=compute_debarcoding(handles);
%% end 

% compute mahalanobis distances
handles=compute_mahal(handles);

handles=compute_well_abundances(handles);

% plot well abundances
set(handles.plottype,'SelectedObject',handles.separation_plot)
handles=separation_plot(handles);

drawnow
set(handles.parent,'pointer','arrow')
drawnow

set(handles.x_popup,'value',1)
set(handles.y_popup,'value',2)

if any(cellfun(@isempty,handles.m(handles.bc_cols)))
    set(handles.x_popup,'string',handles.c(handles.bc_cols))
    set(handles.y_popup,'string',handles.c(handles.bc_cols))
else
    set(handles.x_popup,'string',handles.m(handles.bc_cols))
    set(handles.y_popup,'string',handles.m(handles.bc_cols))
end

handles.leg=cell(1,length(handles.masses));
for i=1:length(handles.masses)
    m_i=handles.m{handles.bc_cols(i)};
    if ~isempty(m_i)
        handles.leg{i}=m_i;
    else
        handles.leg{i}=handles.c{handles.bc_cols(i)};
    end
end


function handles=calculate_cofactors(handles)
% determine a cofactor for each bc channel by pooling the negative barcodes
% for that channel across the populations 

temp_bcind=handles.bcind;
temp_bcind(handles.deltas<0.3)=0;
N=length(temp_bcind);

neg_bcs=cell(1,handles.num_masses);
% pos_bcs=cell(1,handles.num_masses);
bc_list=1:handles.num_codes;
neg_cofactor=zeros(1,handles.num_masses);
% pos_med=zeros(1,handles.num_masses);
for i=1:handles.num_masses
    neg_list=bc_list(handles.key(:,i)==0);
%     pos_list=bc_list(handles.key(:,i)==1);
    neg_cells=false(N,1);
%     pos_cells=false(N,1);
    for j=neg_list
        neg_cells=neg_cells | temp_bcind==j;
    end
%     for j=pos_list
%         pos_cells=pos_cells | temp_bcind==j;
%     end
    neg_bcs{i}=handles.bcs(neg_cells,i); %this was already transformed using default cofactor
%     pos_bcs{i}=handles.bcs(pos_cells,i); %this was already transformed using default cofactor
    neg_cofactor(i)=handles.default_cofactor*sinh(prctile(neg_bcs{i},95)); %untransformed to raw data val
%     pos_med(i)=handles.default_cofactor*sinh(median(pos_bcs{i})); %untransformed to raw data val
end

if any(isnan(neg_cofactor))
    warndlg('Check your barcode key. You may have included a barcode metal that is constant across all occupied samples.')
    neg_cofactor(isnan(neg_cofactor))=5;
end

neg_cofactor(neg_cofactor<5)=5; %5 is default minimum
neg_cofactor(neg_cofactor>100)=100; %5 is default minimum
handles.cofactors=neg_cofactor;

% %find median of pos bcs for each barcode sample
% norm_vals=zeros(size(handles.key));
% for i=1:handles.num_codes
%     bci=handles.bcs(handles.bcind==i,:);
%     pos_masses=handles.key(i,:)==1;
%     norm_vals(i,pos_masses)=median(bci(:,pos_masses));
% end

function handles=recofactor(handles)
%retransform bcs and norm_vals from default cofactoring to variable
%cofactoring

cofactored_bcs=zeros(size(handles.bcs));

for i=1:handles.num_masses
cofactored_bcs(:,i)=asinh(handles.default_cofactor*sinh(handles.bcs(:,i))/handles.cofactors(i));
% norm_vals(:,i)=asinh(handles.default_cofactor*sinh(norm_vals(:,i))/neg_cofactor(i)); %may not use this
handles.cofactored_xt(:,i)=bmtrans(handles.raw_xt,handles.cofactors(i));
end
handles.cofactored_xl=handles.cofactored_xt([1 end],:);

handles.cofactored_bcs=cofactored_bcs;
handles.normbcs=normalize_bcs(cofactored_bcs);


function save_text_Callback(hObject, eventdata, handles)
% hObject    handle to save_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_text as text
%        str2double(get(hObject,'String')) returns contents of save_text as a double


% --- Executes during object creation, after setting all properties.
function save_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in x_popup.
function x_popup_Callback(hObject, eventdata, handles)
% hObject    handle to x_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns x_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from x_popup

plot_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function x_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in y_popup.
function y_popup_Callback(hObject, eventdata, handles)
% hObject    handle to y_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns y_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from y_popup

plot_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function y_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minsep_Callback(hObject, eventdata, handles)
% hObject    handle to minsep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minsep as text
%        str2double(get(hObject,'String')) returns contents of minsep as a double



% --- Executes during object creation, after setting all properties.
function minsep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minsep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxsep_Callback(hObject, eventdata, handles)
% hObject    handle to maxsep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxsep as text
%        str2double(get(hObject,'String')) returns contents of maxsep as a double


% --- Executes during object creation, after setting all properties.
function maxsep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxsep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in range_ypopup.
function range_ypopup_Callback(hObject, eventdata, handles)
% hObject    handle to range_ypopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns range_ypopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from range_ypopup


% --- Executes during object creation, after setting all properties.
function range_ypopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to range_ypopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in range_xpopup.
function range_xpopup_Callback(hObject, eventdata, handles)
% hObject    handle to range_xpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns range_xpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from range_xpopup


% --- Executes during object creation, after setting all properties.
function range_xpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to range_xpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in well_popup.
function well_popup_Callback(hObject, eventdata, handles)
% hObject    handle to well_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns well_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from well_popup

plot_button_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function well_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to well_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles=mahal_cutoff_Callback(hObject, eventdata, handles)
% filter barcode assignments by mahalanobis distance threshold


handles.mahal_cutoff_val=str2double(get(handles.mahal_cutoff,'string'));

within_cutoffs=handles.mahal<handles.mahal_cutoff_val & handles.deltas>handles.sep_cutoff;
for i=1:handles.num_codes
    handles.well_yield(i)=handles.sample_ratio*nnz(handles.bcind==i & within_cutoffs);
end

selectedobj=get(handles.plottype,'SelectedObject');
if ~strcmp(get(selectedobj,'string'),'Separation')
    plot_well_yields(handles)
    plot_button_Callback(hObject, eventdata, handles)
end

guidata(handles.parent,handles)

function plot_well_yields(handles)
% plots current well yields given chosen parameters

ch=get(handles.yield_panel,'children');
delete(ch)
yield_axis=subplot(5,1,[2 3 4],'parent',handles.yield_panel);
bar(yield_axis,handles.well_yield,'facecolor',[0 0.5 0.4])
ylabel('Cell count')
set(yield_axis,'xlim',[0 handles.num_codes+1],'xtick',[])
yl=get(yield_axis,'ylim');
text(1:handles.num_codes,-yl(2)/15*ones(1,handles.num_codes),handles.wellLabels,'rotation',315)

total_yield=round(100*sum(handles.well_yield)/size(handles.x,1));
title(yield_axis,['Barcode Yields with Current Filters: ' num2str(total_yield) '% assigned'],'fontsize',12)



% --- Executes during object creation, after setting all properties.
function mahal_cutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mahal_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pop_back.
function pop_back_Callback(hObject, eventdata, handles)
% hmoves to the previous barcode for plotting

popval=get(handles.well_popup,'value');
if popval>1
    set(handles.well_popup,'value',popval-1);
end
plot_button_Callback(hObject, eventdata, handles)


% --- Executes on button press in pop_fwd.
function pop_fwd_Callback(hObject, eventdata, handles)
% moves to the next barcode for plotting

popval=get(handles.well_popup,'value');
if popval < length(handles.wellLabels);
    set(handles.well_popup,'value',popval+1);
end
plot_button_Callback(hObject, eventdata, handles)


function bc_cols=find_bc_cols_by_mass(colnames, masses)
% find the columns of the fcs file that correspond to those listed in the
% barcode key

num_cols=length(colnames);

bc_cols=zeros(1,length(masses));

if isfloat(masses)
    stringmasses=cell(1,length(masses));
    for j=1:length(masses)
        stringmasses{j}=num2str(masses(j));
    end
    masses=stringmasses;
end

for i=1:num_cols
    for j=1:length(masses)
        if ~isempty(strfind(colnames{i},[masses{j} 'D'])) || ~isempty(strfind(colnames{i},[masses{j} ')']))
            bc_cols(j)=i;
        end
    end
end

if any(bc_cols == 0)
    error('not all barcode channels found')
end

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

function cofactor_Callback(hObject, eventdata, handles)
% not in use: transform data based on user-selected cofactor

old_cofactor=handles.cofactor_val;

handles.cofactor_val=str2double(get(handles.cofactor,'string'));

handles.bcs=asinh(old_cofactor*sinh(handles.bcs)/handles.cofactor_val);

handles.normbcs=normalize_bcs(handles.bcs);
guidata(handles.parent,handles)



% --- Executes during object creation, after setting all properties.
function cofactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cofactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
