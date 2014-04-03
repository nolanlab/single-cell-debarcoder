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

% Last Modified by GUIDE v2.5 18-Oct-2012 13:31:07

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

%setup default initial view and paramters
set(handles.ax,'visible','off')
set(handles.biax_panel,'visible','off')
set(handles.color_panel,'visible','off')
set(handles.plottype,'SelectionChangeFcn',{@plot_changefcn,handles})
set(handles.plottype,'SelectedObject',handles.colorplot)
set(handles.color_panel,'SelectedObject',handles.color_mahal)
set(handles.cutoff_text,'string','30')
set(handles.delta_text,'string','0.1')
handles.sep_cutoff=0.1;
handles.parent=gcf;

%remove unwanted toolbar options
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

%update handles
guidata(hObject, handles)

%load barcode key
bc_button_Callback(hObject, eventdata, handles);


% --- Outputs from this function are returned to the command line.
function varargout = single_cell_debarcode_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function plot_changefcn(hObject,eventdata,handles)
%switches between views for the different plotting options
selectedobj = get(handles.plottype,'SelectedObject');

switch get(selectedobj,'string')
    case 'Color'
        set(handles.biax_panel,'visible','off')
        set(handles.color_panel,'visible','off')
    case 'Single Biaxial'
        set(handles.biax_panel,'visible','on')
        set(handles.color_panel,'visible','on')
    case 'All BC Biaxials'
        set(handles.biax_panel,'visible','off')
        set(handles.color_panel,'visible','on')
end
 
% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selectedobj = get(handles.plottype,'SelectedObject');

wellnum=get(handles.well_popup,'value');
mahal_cutoff=str2double(get(handles.mahal_cutoff,'string'));
sep_cutoff=str2double(get(handles.delta_text,'string'));

if sep_cutoff ~= handles.sep_cutoff %sep_cutoff has been updated, so need to update mahals and bc assignments
    handles.sep_cutoff=sep_cutoff;
    num_codes=length(handles.wellLabels);
    handles.mahal=zeros(size(handles.deltas));
    for i=1:num_codes
        in_bc=handles.gated{i} & (handles.deltas > sep_cutoff);
        bci=handles.bcs(in_bc,:);
        if size(bci,1)>num_codes
            handles.mahal(in_bc)=mahal(bci,bci);
        end
    end
end

thiswell_bin = handles.gated{wellnum} & (handles.mahal<mahal_cutoff) & (handles.deltas > sep_cutoff);

if ~isempty(wellnum)
    
    set(gcf,'pointer','watch')
    drawnow;
    switch get(selectedobj,'string')
        case 'Color'
            
            xl=bmtrans([-20, 10000],5);
            [xt,xtl]=transform_ticks(xl,5);
            
            oldax=get(handles.ax_panel,'children');
            delete(oldax)
            
            ax=zeros(1,2);
            ax(1)=subplot(2,1,1,'parent',handles.ax_panel);
  
            thiswell=handles.bcs(thiswell_bin,:);
            if ~isempty(thiswell)
                s=1:size(thiswell,1);
                [ax,h1,h2]=plotyy(s,thiswell,s,thiswell);
                hold on
                
                set(ax,'Box','off','xtick',[],'Layer','top')
                set(h2,'visible','off')
                set(h1,'marker','.','linestyle','none')
                
                set(ax(1),'ytick',xt,'yticklabel',xtl)
                set(get(ax(1),'Ylabel'),'String','untransformed values','fontsize',12)
                set(get(ax(2),'Ylabel'),'String','asinh-transformed values','fontsize',12)
                
                legend(handles.leg,'location','northeastoutside')
            end
            title(num2str(handles.key(wellnum,:)),'fontweight','bold','fontsize',12);
            
            %norm bcs
            ax(2)=subplot(2,1,2,'parent',handles.ax_panel);
            
            thiswell=handles.normbcs(thiswell_bin,:);
            
            if ~isempty(thiswell)
                s=1:size(thiswell,1);
                plot(s,thiswell,'.');
                
                set(get(ax(2),'Ylabel'),'String','normalized values','fontsize',12)
                
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
            
            thiswell=handles.bcs(thiswell_bin,[xcol ycol]);
            mdists=handles.mahal(thiswell_bin);
            seps=handles.deltas(thiswell_bin);
            if ~isempty(thiswell)

                if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                    scatter(thiswell(:,1),thiswell(:,2),4,mdists)
                    set(gca,'clim',[0 mahal_cutoff])
                    colormap(hsv)
                else
                    scatter(thiswell(:,1),thiswell(:,2),4,seps)
                    set(gca,'clim',[sep_cutoff 0.5])
                    colormap(flipud(hsv))
                end
                
                cb=colorbar;
                if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                    set(get(cb,'ylabel'),'string','Mahalanobis Distance','fontsize',14)
                else
                    set(get(cb,'ylabel'),'string','Separation','fontsize',14)
                end
                
            end
            xl=bmtrans([-20, 10000],5);
            [xt,xtl]=transform_ticks(xl,5);
            set(gca,'xlim',xl,'ylim',xl,'xtick',xt,'xticklabel',xtl,'ytick',xt,'yticklabel',xtl)
            title(num2str(handles.key(wellnum,:)),'fontweight','bold','fontsize',12);
            
        case 'All BC Biaxials'
            oldax=get(handles.ax_panel,'children');
            delete(oldax)
            n=size(handles.bcs,2);
            xl=bmtrans([-20, 10000],5);
            [xt,xtl]=transform_ticks(xl,5);
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
                        
                        thiswell=handles.bcs(thiswell_bin,[j i+1]);
                        mdists=handles.mahal(thiswell_bin);
                        seps=handles.deltas(thiswell_bin);
                        
                        if ~isempty(thiswell)
                            if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                                scatter(thiswell(:,1),thiswell(:,2),2,mdists)
                                set(gca,'clim',[0 mahal_cutoff])
                                colormap(hsv)
                            else
                                scatter(thiswell(:,1),thiswell(:,2),2,seps)
                                set(gca,'clim',[sep_cutoff 0.5])
                                colormap(flipud(hsv))
                            end
                        end                        
                        set(ax,'xlim',xl,'ylim',xl,'xtick',[],'ytick',[])
                    elseif j>0 && i==j-1 %&& any([i j])
                        ax=subplot(n+1,n+1,Ind,'Parent',handles.ax_panel);
                        if sum(handles.gated{wellnum}) ~= 0
                            [binsize,binloc]=hist(handles.bcs(handles.gated{wellnum} & handles.mahal<mahal_cutoff,j),100);
                            bar(binloc,binsize,'edgecolor','none','facecolor',[0 0.5 0])
                        end
                        set(ax,'xlim',xl,'xtick',[],'ytick',[])                       
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
                set(bigax,'clim',[0 mahal_cutoff])
                cb=colorbar('position',[0.85 0.5 0.027 0.4]);
                set(get(cb,'ylabel'),'string','Mahalanobis Distance','fontsize',14)
            else
                set(bigax,'clim',[sep_cutoff 0.5])
                cb=colorbar('position',[0.85 0.5 0.027 0.4]);
                set(get(cb,'ylabel'),'string','Separation','fontsize',14)
            end
            set(bigax,'visible','off')
    end
    
    drawnow
    set(gcf,'pointer','arrow')
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


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


FileName=get(handles.save_text,'string');

mahal_cutoff=str2double(get(handles.mahal_cutoff,'string'));
sep_cutoff=str2double(get(handles.delta_text,'string'));


PathName = uigetdir;
if PathName ~= 0
    
    set(gcf,'pointer','watch')
    drawnow
    
    fid=fopen([PathName '/' FileName '_BarCodeList.txt'],'w');  %this is to record the barcode labels for use later with spade
    
    wellNumbers=zeros(size(handles.x,1),1); %%%
    m=[handles.m 'barcode']; %add barcode column
    c=[handles.c 'barcode'];
    n=length(c);
    
    %% need to recompute handles.gated using all (not just sampled) cells
    %%Aug2013
    handles.bcs=bmtrans(handles.x(:,handles.bc_cols),5); %switching sampled bcs to full bcs
    
    percs=prctile(handles.bcs,[1 99]);
    ranges=diff(percs);
    deltas=bsxfun(@minus,handles.bcs,percs(1,:));
    handles.normbcs=bsxfun(@rdivide,deltas,ranges);
    
    if length(unique(sum(handles.key,2)))==1
        
        [sorted,ix]=sort(handles.normbcs,2,'descend');
        
        % to look at top k barcodes, rather than barcodes above largest separation
        numdf=sum(handles.key(1,:));  %this is how many positive we are expecting
        
        lowests=sorted(:,numdf);
        
        %need to get rid of ones that aren't above actual signal cutoff
        N=size(handles.bcs,1);
        indlist=(1:N)';
        inds1=sub2ind(size(ix),indlist,ix(:,numdf));
        cutoff1=bmtrans(str2double(get(handles.cutoff_text,'string')),5);
        toolow1=handles.bcs(inds1)<cutoff1;
        
        lowests(toolow1)=nan;
        
        
        deltas=sorted(:,numdf)-sorted(:,numdf+1);  %separation between 3rd and 4th highest normalized barcodes
        
    else
        
        [sorted,ix]=sort(handles.normbcs,2,'ascend');
        
        seps=diff(sorted,1,2);
        
        N=size(handles.bcs,1);
        indlist=(1:N)';
        
        [~,locs]=sort(seps,2,'descend');  %locs are locations in ix of bc level that is on lower side of largest gap.  i.e., if locs is 5, the largest bc separation is between barcode ix(5) and ix(6)
        
        cutoff=bmtrans(str2double(get(handles.cutoff_text,'string')),5);
        
        indsabove=sub2ind(size(ix),indlist,locs(:,1)+1);  %+1 so that finding the lowest "real" bc
        betws=ix(indsabove);
        indsabove=sub2ind(size(handles.bcs),indlist,betws);
        % lowests=handles.bcs(inds);  %transformed values of lowest "good" bcs  %switched to normalized bcs
        lowests=handles.normbcs(indsabove);  %normalized transformed values of lowest "good" bcs
        
        %%%finding the highest "non-real" bc
        indsbelow=sub2ind(size(ix),indlist,locs(:,1));
        betws=ix(indsbelow);
        indsbelow=sub2ind(size(handles.bcs),indlist,betws);
        % nextlowests=handles.bcs(inds); %switched to normalized bcs
        nextlowests=handles.normbcs(indsbelow);
        
        %%%
        
        toolow=find(handles.bcs(indsabove)<cutoff);  %these aren't high enough to count.  go to next-biggest-sep. using actual bcs here, not normalized
        inds=sub2ind(size(ix),toolow,locs(toolow,2)+1);
        betws2=ix(inds);
        inds=sub2ind(size(handles.bcs),toolow,betws2);
        lowests2=handles.normbcs(inds);
        highernow=handles.bcs(inds)>cutoff;  %again using actualy bcs, not normalized, to check against cutoff
        %might still need to account for when the largest sep is high ...  can
        %first try eliminating these just with illegal barcodes
        lowests(toolow(highernow))=lowests2(highernow);
        lowests(toolow(~highernow))=nan;  %if the second try didn't find one >10, set to nan
        
        %adding in the replaced "non-real" bc
        inds=sub2ind(size(ix),toolow,locs(toolow,2));
        betws2=ix(inds);
        inds=sub2ind(size(handles.bcs),toolow,betws2);
        % nextlowests2=handles.bcs(inds);  %switched to normalized bcs
        nextlowests2=handles.normbcs(inds);
        nextlowests(toolow(highernow))=nextlowests2(highernow);
        nextlowests(toolow(~highernow))=nan;
        
        deltas=lowests-nextlowests;
        
    end
    
    handles.deltas=deltas;
    
    code_assign=false(N,length(handles.masses));  %rows will be barcodes of each cell
    for j=1:length(handles.masses)
        code_assign(:,j)=handles.normbcs(:,j) >= lowests;
    end
    
    num_codes=length(handles.wellLabels);
    length_codes=size(handles.key,2);
    handles.gated=cell(1,num_codes);
    num_cells=size(handles.bcs,1);
    
    
    for i=1:num_codes
        clust_inds=true(num_cells,1);
        for j=1:length_codes
            clust_inds=clust_inds & (code_assign(:,j)==handles.key(i,j));
        end
        handles.gated{i}=clust_inds;
    end
    
    
    %%% temporary to look at trimming
    
    handles.mahal=zeros(size(handles.deltas));
    for i=1:num_codes
        in_bc=handles.gated{i} & (handles.deltas > sep_cutoff);
        bci=handles.bcs(in_bc,:);
        if size(bci,1)>num_codes
            handles.mahal(in_bc)=mahal(bci,bci);
        end
    end
    
    %% end re-computation of debarcoding
    
    not_inawell=true(size(handles.gated{1}));
    
    
    for i=1:length(handles.wellLabels)
        
        fprintf(fid,'%s\n',handles.wellLabels{i});  %printing out one line of _BarCodeList.txt
        
        thiswell_bin = handles.gated{i} & (handles.mahal<mahal_cutoff) & (handles.deltas > sep_cutoff);
        not_inawell(thiswell_bin)=false; %cells is this well removed from unassigned_binary
        
        data=zeros(sum(thiswell_bin),n);
        data(:,1:end-1)=handles.x(thiswell_bin,:);
        data(:,end)=i;
        if ~isempty(data)
            fca_writefcs([PathName '/' FileName '_' handles.wellLabels{i} '.fcs'],data,m,c)
            wellNumbers(thiswell_bin)=i; %%%
        end
    end
    
    fclose(fid);
    
    data=zeros(sum(not_inawell),n);
    data(:,1:end-1)=handles.x(not_inawell,:);
    
    fca_writefcs([PathName '/' FileName '_unassigned.fcs'],data,m,c)
    %     fca_writefcs([PathName '/' FileName '_BClabeled.fcs'],[handles.x wellNumbers],[handles.m 'barcode'],[handles.c 'barcode']) %%%
    %     dlmwrite([PathName '/' FileName '_wellNumbers.txt'],wellNumbers,'delimiter','\n')
    
    fid=fopen([PathName '/Debarcode_Parameters.txt'],'w');
    fprintf(fid,'%s\n',datestr(now));
    fprintf(fid,'%s\n','Input files:');
    for i=1:length(handles.current_files)
        fprintf(fid,'\t%s\n',handles.current_files{i});
    end
    fprintf(fid,'%s\n',['Barcode Key: ' handles.key_filename]);
    fprintf(fid,'%s\n',['Mahalanobis cutoff: ' get(handles.mahal_cutoff,'string')]);
    fprintf(fid,'%s\n',['Separation cutoff: ' get(handles.delta_text,'string')]);
    fprintf(fid,'%s\n',['Output saved to: ' PathName '/']);
    fclose(fid);
    
    set(gcf,'pointer','arrow')
    drawnow
else
    return
end


% --- Executes on button press in folder_button.
function folder_button_Callback(hObject, eventdata, handles)
% hObject    handle to folder_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[files,handles.pathname] = uigetfile('*.fcs','Select fcs files to debarcode','multiselect','on');

if handles.pathname ~= 0
    
    set(gcf,'pointer','watch')
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
    
    guidata(handles.parent,handles)
    
    load_bc_data(handles)
    
    ch=get(handles.ax_panel,'children');
    delete(ch)
    drawnow
    set(gcf,'pointer','arrow')
    
else
    return
end



% --- Executes on button press in debarcode_button.
function debarcode_button_Callback(hObject, eventdata, handles)
% hObject    handle to debarcode_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'x')
    msgbox('You must first load fcs files!')
    return
end

if ~isfield(handles,'bcs')
    msgbox('Check your barcode key.')
    return
end

ch=get(handles.ax_panel,'children');
delete(ch)
handles.ax=axes('parent',handles.ax_panel);
set(gcf,'currentaxes',handles.ax)
set(handles.ax,'visible','off')
tb=text(0.5,0.5,'Debarcoding ...','horizontalalignment','center','verticalalignment','middle','fontsize',14);

set(gcf,'pointer','watch')
drawnow

% [sorted,ix]=sort(handles.bcs,2); %switched to normalized bcs



if length(unique(sum(handles.key,2)))==1
    
    [sorted,ix]=sort(handles.normbcs,2,'descend');
    
    % to look at top k barcodes, rather than barcodes above largest separation
    numdf=sum(handles.key(1,:));  %this is how many positive we are expecting
    
    
    lowests=sorted(:,numdf);
    
    %need to get rid of ones that aren't above actual signal cutoff
    N=size(handles.bcs,1);
    indlist=(1:N)';
    inds1=sub2ind(size(ix),indlist,ix(:,numdf));
    cutoff1=bmtrans(str2double(get(handles.cutoff_text,'string')),5);
    toolow1=handles.bcs(inds1)<cutoff1;
    
    lowests(toolow1)=nan;
    
    deltas=sorted(:,numdf)-sorted(:,numdf+1);  %separation between 3rd and 4th highest normalized barcodes
    
else
    
    [sorted,ix]=sort(handles.normbcs,2,'ascend');
    
    seps=diff(sorted,1,2);
    
    N=size(handles.bcs,1);
    indlist=(1:N)';
    
    [~,locs]=sort(seps,2,'descend');  %locs are locations in ix of bc level that is on lower side of largest gap.  i.e., if locs is 5, the largest bc separation is between barcode ix(5) and ix(6)
    
    cutoff=bmtrans(str2double(get(handles.cutoff_text,'string')),5);
    
    indsabove=sub2ind(size(ix),indlist,locs(:,1)+1);  %+1 so that finding the lowest "real" bc
    betws=ix(indsabove);
    indsabove=sub2ind(size(handles.bcs),indlist,betws);
    % lowests=handles.bcs(inds);  %transformed values of lowest "good" bcs  %switched to normalized bcs
    lowests=handles.normbcs(indsabove);  %normalized transformed values of lowest "good" bcs
    
    %%%finding the highest "non-real" bc
    indsbelow=sub2ind(size(ix),indlist,locs(:,1));
    betws=ix(indsbelow);
    indsbelow=sub2ind(size(handles.bcs),indlist,betws);
    % nextlowests=handles.bcs(inds); %switched to normalized bcs
    nextlowests=handles.normbcs(indsbelow);
    
    %%%
    
    toolow=find(handles.bcs(indsabove)<cutoff);  %these aren't high enough to count.  go to next-biggest-sep. using actual bcs here, not normalized
    inds=sub2ind(size(ix),toolow,locs(toolow,2)+1);
    betws2=ix(inds);
    inds=sub2ind(size(handles.bcs),toolow,betws2);
    % lowests2=handles.bcs(inds); %switched to normalized bcs
    lowests2=handles.normbcs(inds);
    highernow=handles.bcs(inds)>cutoff;  %again using actualy bcs, not normalized, to check against cutoff
    %might still need to account for when the largest sep is high ...  can
    %first try eliminating these just with illegal barcodes
    lowests(toolow(highernow))=lowests2(highernow);
    lowests(toolow(~highernow))=nan;  %if the second try didn't find one >10, set to nan
    
    %adding in the replaced "non-real" bc
    inds=sub2ind(size(ix),toolow,locs(toolow,2));
    betws2=ix(inds);
    inds=sub2ind(size(handles.bcs),toolow,betws2);
    % nextlowests2=handles.bcs(inds);  %switched to normalized bcs
    nextlowests2=handles.normbcs(inds);
    nextlowests(toolow(highernow))=nextlowests2(highernow);
    nextlowests(toolow(~highernow))=nan;
    
    deltas=lowests-nextlowests;
    
end


%don't get rid of anything based on separatio now; save sep value
%and use like mahal value

%         lowests(deltas<str2double(get(handles.delta_text,'string')))=nan;  %if the largest separation is less than sep cutoff in asinh space, dont count it
handles.deltas=deltas;

code_assign=false(N,length(handles.masses));  %rows will be barcodes of each cell
for j=1:length(handles.masses)
    % code_assign(:,j)=handles.bcs(:,j) >= lowests; %switched to normalized bcs
    code_assign(:,j)=handles.normbcs(:,j) >= lowests; %this takes the min cutoff into account
end

% gated=gate_bcs(code_assign,handles.x); %%replace with below using csv key
num_codes=length(handles.wellLabels);
length_codes=size(handles.key,2);
handles.gated=cell(1,num_codes);
num_cells=size(handles.bcs,1);

%         handles.dump=true(num_cells,1);

for i=1:num_codes
    clust_inds=true(num_cells,1);
    for j=1:length_codes
        clust_inds=clust_inds & (code_assign(:,j)==handles.key(i,j));
    end
    handles.gated{i}=clust_inds;
    %             handles.dump(clust_inds)=false;  %these cells are in a gate
end

%%% temporary to look at trimming
handles.mahal=zeros(size(handles.deltas));
for i=1:num_codes    
    bci=handles.bcs(handles.gated{i},:);
    if size(bci,1)>num_codes
        handles.mahal(handles.gated{i})=mahal(bci,bci);
    end
end

%%% plot well abundances

numseps=10;
minsep=0;
maxsep=0.5;
handles.seprange=linspace(minsep,maxsep,numseps);

clust_size=zeros(numseps,num_codes);

for i=1:numseps
    for j=1:length(handles.wellLabels)
        clust_size(i,j) = sum(handles.gated{j} & (handles.deltas > handles.seprange(i)));
    end
end

delete(tb)

handles.ax=axes('parent',handles.ax_panel);
handles.lines=plot(handles.ax,handles.seprange,clust_size);
set(get(handles.ax,'XLabel'),'String','Minimum separation','fontsize',12)
set(get(handles.ax,'YLabel'),'String','Number of cells','fontsize',12)
ylabel('Number of cells','fontsize',12)

set(handles.lines,'ButtonDownFcn',{@select_line,handles});

drawnow
set(gcf,'pointer','arrow')
drawnow

%         mahal_threshold=30;
%
%             handles.bclabels=zeros(size(handles.bcs,1),1);
%             for i=1:num_codes
%                 handles.bclabels(handles.gated{i})=i;
%             end
%             inwell=find(handles.bclabels~=0);
%             bclabels=handles.bclabels(inwell);
%             bcdata=handles.bcs(inwell,:);
%             gmfit=gmdistribution.fit(bcdata, num_codes,'start',bclabels);
%             mdist=mahal(gmfit,bcdata);
% %             totrim=false(size(inwell));
% %             for i=1:num_codes
% %                 totrim(bclabels == i & mdist(:,i) > mahal_threshold) = true;
% %             end
% %             handles.bclabels(inwell(totrim))=0;
%             handles.mahal=zeros(size(handles.bcs,1),num_codes);
%             handles.mahal(inwell,:)=mdist;
%          %%%




guidata(hObject,handles)


function select_line(button,eventdata,handles)

ch=get(handles.ax,'children');
T=ch(strcmp('text',get(ch,'type')));
if ~isempty(T)
    delete(T)
end

ind=handles.lines==gco;
cp=get(handles.ax,'currentpoint');
text(cp(1,1),cp(1,2),handles.wellLabels{ind},'fontsize',12,'verticalalignment','baseline','edgecolor','k','backgroundcolor','w','interpreter','none')



function delta_text_Callback(hObject, eventdata, handles)
% hObject    handle to delta_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delta_text as text
%        str2double(get(hObject,'String')) returns contents of delta_text as a double


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
% hObject    handle to bc_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.key_filename, handles.bcpathname]=uigetfile({'*.csv','*.CSV'},'Select barcode key');

if handles.bcpathname ~= 0
    
    x=importdata([handles.bcpathname handles.key_filename]);
    handles.masses=cellstr(num2str(x.data(1,:)'));
    str=sprintf('Using key:\n%s',handles.key_filename);
    set(handles.key_text,'string',str);
    
    handles.masses=cellstr(num2str(x.data(1,:)'));
    handles.wellLabels=x.textdata(2:end);
    set(handles.well_popup,'string',handles.wellLabels)
    
    handles.key=x.data(2:end,:);
    
    guidata(handles.parent,handles)
    
    %if fcs file already loaded, update which are the bc cols and the bc data i
    if isfield(handles,'c')
        load_bc_data(handles)
    end
    
else
    return
end

function load_bc_data(handles)

try
    handles.bc_cols=find_bc_cols_by_mass(handles.c,handles.masses);
catch err
    if strcmp('not all barcode channels found',err.message)
        set(gcf,'pointer','arrow')
        msgbox('Check your barcode key.')
        guidata(handles.parent,handles)
        return
    else
        set(gcf,'pointer','arrow')
        rethrow(err)
    end
end

%sample 100000 cells and use until save
num_cells=size(handles.x,1);
sample_size=100000;
if num_cells>sample_size
    handles.bcs=bmtrans(handles.x(randsample(num_cells,sample_size),handles.bc_cols),5);
else
    handles.bcs=bmtrans(handles.x(:,handles.bc_cols),5);  %matrix of each cell's bc channels, transformed
end
%

percs=prctile(handles.bcs,[1 99]);
ranges=diff(percs);  %difference between 99th and 1st percentile of bc channels
deltas=bsxfun(@minus,handles.bcs,percs(1,:));
handles.normbcs=bsxfun(@rdivide,deltas,ranges);  %could still collect on edges (0,1)

set(handles.x_popup,'value',1)
set(handles.y_popup,'value',1)

set(handles.x_popup,'string',handles.m(handles.bc_cols))
set(handles.y_popup,'string',handles.m(handles.bc_cols))

handles.leg=cell(1,length(handles.masses));
for i=1:length(handles.masses)
    handles.leg{i}=handles.m{handles.bc_cols(i)};
end

guidata(handles.parent,handles)



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



function mahal_cutoff_Callback(hObject, eventdata, handles)
% hObject    handle to mahal_cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mahal_cutoff as text
%        str2double(get(hObject,'String')) returns contents of mahal_cutoff as a double


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
% hObject    handle to pop_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

popval=get(handles.well_popup,'value');
if popval>1
    set(handles.well_popup,'value',popval-1);
end
plot_button_Callback(hObject, eventdata, handles)


% --- Executes on button press in pop_fwd.
function pop_fwd_Callback(hObject, eventdata, handles)
% hObject    handle to pop_fwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

popval=get(handles.well_popup,'value');
if popval < length(handles.wellLabels);
    set(handles.well_popup,'value',popval+1);
end
plot_button_Callback(hObject, eventdata, handles)

function bc_cols=find_bc_cols_by_mass(colnames, masses)

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

function [xt,xtl]=transform_ticks(xL,cofactor)

xticks=[];
xtl=[];

if xL(1)<0
    xLa=ceil(log10(abs(cofactor*sinh(xL(1)))));
    xLb=ceil(log10(cofactor*sinh(xL(2))));
    ints=-1:0.1:-0.2;
    
    Ind=xLa;
    for i=1:xLa
        xticks=[xticks ints*10^Ind];
        xtl=[xtl; num2str(-10^Ind); cell(8,1)];
        Ind=Ind-1;
    end
    xticks=[xticks -1 0 1];
    xtl=[xtl; cell(1,1); '0'; cell(1,1)];
    
    ints=0.2:0.1:1;
    
    for i=1:xLb
        xticks=[xticks ints*10^i];
        xtl=[xtl; cell(8,1); num2str(10^i)];
    end
else
    xLb=ceil(log10(cofactor*sinh(xL(2))));
    xticks=[xticks  0 1];
    xtl=['0'; cell(1,1)];
    
    ints=0.2:0.1:1;
    
    for i=1:xLb
        xticks=[xticks ints*10^i];
        xtl=[xtl; cell(8,1); num2str(10^i)];
    end
    
end

xt=bmtrans(xticks,cofactor);

function denslevel=densassign(x,perc,M)
%denslevel=densassign(x,perc,M)

if nargin<3
    M=128;
end

if nargin<2
    perc=0;
end

n=size(x,1);  %number of data points
d=size(x,2);  %number of dimensions
MM = M^d; %number of total gridpoints

mins=min(x);  %vector of min of each column of data
maxs=max(x);   %vector of max of each column of data
Diff=maxs-mins;
mins=mins-perc/100*Diff;
maxs=maxs+perc/100*Diff;

Delta = 1/(M-1)*(maxs-mins);  %vector of distances between neighbor grid points in each dimension

% ym=gridpoints(d,1:M); %list of every possible combo of 1:M

ye=zeros(d,M);
multby=zeros(1,d);  %used in coord transfrom from m to k
pointLL=zeros(n,d);  %this will be the "lower left" gridpoint to each data point
for i = 1:d
    ye(i,:) = linspace(mins(i),maxs(i),M);
    multby(i)=M^(i-1);
    pointLL(:,i)=floor((x(:,i)-mins(i))./Delta(i)) + 1;
end
pointLL(pointLL==M)=M-1;  %this avoids going over grid boundary

%% assign each data point to its closest grid point
[xgrid,ygrid]=meshgrid(ye(1,:),ye(2,:));
z=reshape(1:MM,M,M);
GridAssign=interp2(xgrid,ygrid,z',x(:,1),x(:,2),'nearest');  %this associates each data point with its nearest grid point

%% compute w
Deltmat=repmat(Delta,n,1);
shape=M*ones(1,d);

wmat=zeros(M,M);
for i=0:1  %number of neighboring gridpoints in d dimensions
    for j=0:1
        pointm=pointLL+repmat([j i],n,1);  %indices of ith neighboring gridpoints
        pointy=zeros(n,d);
        for k=1:d
            pointy(:,k)=ye(k,pointm(:,k));  %y-values of ith neighboring gridpoints
        end
        W=prod(1-(abs(x-pointy)./Deltmat),2);  %contribution to w from ith neighboring gridpoint from each datapoint
        wmat=wmat+accumarray(pointm,W,shape);  %sums contributions for ith gridpoint over data points and adds to wmat
    end
end

%% compute f, sig, df and A
n6 = n^(-1/6);
h=zeros(1,d);
Z=zeros(1,d);
Zin=cell(1,d);
for i =1:d
    h(i) = std(x(:,i))*n6;
    Z(i) = min(floor(4*h(i)/Delta(i)), M-1);
    Zin{i}=-Z(i):Z(i);
end

phi = @(x) 1/sqrt(2*pi)*exp(-x.^2./2);

[L{1},L{2}]=meshgrid(Zin{1},Zin{2});

Phix=phi(L{1}*Delta(1)./h(1))./h(1);
Phiy=phi(L{2}*Delta(2)./h(2))./h(2);

Phimat = (Phix.*Phiy)';   %matrix of Phi for inputting into convn

fmat = 1/n*conv2(wmat,Phimat,'same');  %d-dim matrix of estimated densities
f=reshape(fmat,1,MM);  %f in single vector

[~,bin]=histc(f,linspace(0,max(f),64));
denslevel=bin(GridAssign);

function y=bmtrans(x,c)
%c=5;
y=asinh(1/c*x);
