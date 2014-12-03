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
% handles.obj.sep_cutoff=0.3; %now in obj

set(handles.mahal_cutoff,'string','30')
% handles.obj.mahal_cutoff_val=30; %now in obj

% handles.obj.default_cofactor=5; %now in obj

% axis limits and ticks %now in obj
% axticks=load('axticks.mat');
% handles.obj.xtl=axticks.xtl;
% handles.obj.raw_xl=[-10, 10000];
% handles.obj.raw_xt=axticks.xt;
% handles.obj.default_xl=bmtrans(handles.obj.raw_xl,handles.obj.default_cofactor);
% handles.obj.default_xt=bmtrans(axticks.xt,handles.obj.default_cofactor);

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
thiswell_bin = (handles.obj.bcind==wellnum) & (handles.obj.mahal<handles.obj.mahal_cutoff_val) & (handles.obj.deltas > handles.obj.sep_cutoff);

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
            
            thiswell=handles.obj.bcs(thiswell_bin,:);
            if ~isempty(thiswell)
                s=1:size(thiswell,1);
                plot3(s,thiswell,rand(size(s))-1,'.')
                hold on              
                set(ax(1),'Box','off','xtick',[],'Layer','top')
                set(ax(1),'ytick',handles.obj.default_xt,'yticklabel',handles.obj.xtl)
                set(get(ax(1),'Ylabel'),'String','untransformed values','fontsize',12)
                set(get(ax(1),'Ylabel'),'String','Barcode intensities','fontsize',12)
                xlabel('Event','fontsize',12)
                legend(handles.leg,'location','northeastoutside')    
            end
            title(num2str(handles.obj.key(wellnum,:)),'fontweight','bold','fontsize',12);
            
            %norm bcs
            ax(2)=subplot(2,1,2,'parent',handles.ax_panel);
            hold on
            
            thiswell=handles.obj.normbcs(thiswell_bin,:);           
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
            
            thiswell=handles.obj.cofactored_bcs(thiswell_bin,[xcol ycol]);
            
            mdists=handles.obj.mahal(thiswell_bin);
            
            seps=handles.obj.deltas(thiswell_bin);
            if ~isempty(thiswell)
                
                if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                    scatter(thiswell(:,1),thiswell(:,2),4,mdists)
                    set(gca,'clim',[0 handles.obj.mahal_cutoff_val])
                    colormap(flipud(cm))
                else
                    scatter(thiswell(:,1),thiswell(:,2),4,seps)
                    set(gca,'clim',[handles.obj.sep_cutoff 1])
                    
                    colormap(cm)
                end
                
                cb=colorbar;
                if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                    set(get(cb,'ylabel'),'string','Mahalanobis Distance','fontsize',14)
                else
                    set(get(cb,'ylabel'),'string','Separation','fontsize',14)
                end
                
            end
            
            set(gca,'xlim',handles.obj.cofactored_xl(:,xcol),'ylim',handles.obj.cofactored_xl(:,ycol),...
                'xtick',handles.obj.cofactored_xt(:,xcol),'xticklabel',handles.obj.xtl,...
                'ytick',handles.obj.cofactored_xt(:,ycol),'yticklabel',handles.obj.xtl)
            
            title(num2str(handles.obj.key(wellnum,:)),'fontweight','bold','fontsize',12);
            
        case 'All Barcode Biaxials'
            oldax=get(handles.ax_panel,'children');
            delete(oldax)
            n=size(handles.obj.bcs,2);
            
            cm=jet(64);
            
            Ind=1;
            for i=0:n
                for j=0:n
                    
                    if j==0 && i~=n
                        ax=subplot(n+1,n+1,Ind,'Parent',handles.ax_panel);
                        set(ax,'visible','off')
                        text(0.5,0.5,handles.obj.masses{i+1},'horizontalalignment','center','verticalalignment','middle')
                    elseif j<=i && i<n
                        ax=subplot(n+1,n+1,Ind,'Parent',handles.ax_panel);
                        hold on

                        thiswell=handles.obj.cofactored_bcs(thiswell_bin,[j i+1]);

                        mdists=handles.obj.mahal(thiswell_bin);
                        seps=handles.obj.deltas(thiswell_bin);
                        
                        if ~isempty(thiswell)
                            if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                                scatter(thiswell(:,1),thiswell(:,2),2,mdists)
                                set(gca,'clim',[0 handles.obj.mahal_cutoff_val])
                            else
                                scatter(thiswell(:,1),thiswell(:,2),2,seps)
                                set(gca,'clim',[handles.obj.sep_cutoff 1])
                            end
                        end
                        
                        set(ax,'xlim',handles.obj.cofactored_xl(:,j),'ylim',handles.obj.cofactored_xl(:,i+1),'xtick',[],'ytick',[])
                    elseif j>0 && i==j-1 
                        ax=subplot(n+1,n+1,Ind,'Parent',handles.ax_panel);
                        if nnz(handles.obj.bcind==wellnum) ~= 0
                            [binsize,binloc]=hist(handles.obj.cofactored_bcs(thiswell_bin,j),100);
                            bar(binloc,binsize,'edgecolor','none','facecolor',[0 0.5 0])
                        end
                        set(ax,'xlim',handles.obj.cofactored_xl(:,j),'xtick',[],'ytick',[])
                    elseif i==n && j>0
                        ax=subplot(n+1,n+1,Ind,'Parent',handles.ax_panel);
                        set(ax,'visible','off')
                        text(0.5,0.5,handles.obj.masses{j},'horizontalalignment','center','verticalalignment','middle')
                    end
                    Ind=Ind+1;
                end
            end
            
            bigax=axes('parent',handles.ax_panel);
            if get(handles.color_panel,'SelectedObject') == handles.color_mahal
                colormap(flipud(cm))
                set(bigax,'clim',[0 handles.obj.mahal_cutoff_val])
                cb=colorbar('position',[0.85 0.5 0.027 0.4]);
                set(get(cb,'ylabel'),'string','Mahalanobis Distance','fontsize',14)
            else
                colormap(cm)
                set(bigax,'clim',[handles.obj.sep_cutoff 1])
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

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


FileName=get(handles.save_text,'string');

PathName = uigetdir;

if PathName ~= 0
    
    set(handles.parent,'pointer','watch')
    drawnow
    
    handles.obj.write_bc_fcs_files(PathName,FileName)

    % write log file of debarcoding process
    fid=fopen([PathName filesep 'Debarcode_Parameters.txt'],'w');
    fprintf(fid,'%s\n',datestr(now));
    fprintf(fid,'%s\n','Input files:');
    for i=1:length(handles.current_files)
        fprintf(fid,'\t%s\n',handles.current_files{i});
    end
    fprintf(fid,'%s\n',['Barcode Key: ' handles.obj.key_filename]);
    fprintf(fid,'%s\n',['Mahalanobis cutoff: ' get(handles.mahal_cutoff,'string')]);
    fprintf(fid,'%s\n',['Separation cutoff: ' get(handles.delta_text,'string')]);
    fprintf(fid,'%s\n',['Output saved to: ' PathName filesep]);
    fclose(fid);
    
    set(handles.parent,'pointer','arrow')
    drawnow
else
    return
end

% --- Executes on button press in folder_button.
function folder_button_Callback(hObject, eventdata, handles)
% read in an fcs file. if multiple files are selected, they are
% concatenated

[files,pathname] = uigetfile('*.fcs','Select fcs files to debarcode','multiselect','on');

if pathname ~= 0 %didn't hit cancel   
    set(handles.parent,'pointer','watch')
    drawnow
    
    if iscell(files) %more than 1 file selected, which means concatenation
        filenames=strcat(pathname,files);
        formatstr=repmat('\n%s',[1 length(files)]);
        str=sprintf(['Using files:' formatstr],files{:});
        handles.current_files=files;
    else
        filenames=fullfile(pathname,files);
        str=sprintf('Using file:\n%s',files);
        handles.current_files={files};
    end
    set(handles.file_text,'string',str);
    
    handles.obj=handles.obj.load_fcs_files(filenames);
    
    % find bc columns and do preliminary transformation
    handles=load_bc_data(hObject, eventdata,handles);
    
    guidata(hObject,handles)
    
    set(handles.parent,'pointer','arrow')
    drawnow
else
    return
end

function handles=separation_plot(handles)
% makes a histogram of total yield binned by barcode separation, and a plot of each well's yield as
% a function of separation cutoff

ch=get(handles.yield_panel,'children');
delete(ch)
yield_axis=axes('parent',handles.yield_panel);
[hi,xi]=hist(handles.obj.deltas,100);
bar(yield_axis,xi,handles.obj.sample_ratio*hi,'facecolor',[0 0.5 0.4])
set(get(yield_axis,'XLabel'),'String','Barcode separation','fontsize',12)
set(get(yield_axis,'YLabel'),'String','Event count','fontsize',12)
set(yield_axis,'xlim',[0 1],'box','off')

ch=get(handles.ax_panel,'children');
delete(ch)
handles.ax=axes('parent',handles.ax_panel);
set(handles.ax,'colororder',flipud(jet(20)))
hold on
handles.lines=plot(handles.ax,handles.obj.seprange,handles.obj.sample_ratio*handles.obj.clust_size);
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
text(cp(1,1),cp(1,2),handles.obj.wellLabels{ind},'fontsize',12,'verticalalignment','baseline','edgecolor','k','backgroundcolor','w','interpreter','none')


function handles=delta_text_Callback(hObject, eventdata, handles)
% executes upon a change to the separation cutoff value

% update separation cutoff in scd object whenever it is changed in GUI
handles.obj.sep_cutoff=str2double(get(handles.delta_text,'string'));

% compute mahalanobis distances
handles.obj=handles.obj.compute_mahal;

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

[key_filename, key_path]=uigetfile({'*.csv','*.CSV'},'Select barcode key');

handles.obj=scd(fullfile(key_path,key_filename));

if ~isempty(handles.obj.key_filename)

   str=sprintf('Using key:\n%s',handles.obj.key_filename);
   set(handles.key_text,'string',str);
  
    %if fcs file already loaded, update which are the bc cols and the bc data i
%     if isfield(handles,'c')
%         handles=load_bc_data(hObject, eventdata,handles);
%     end
    guidata(hObject,handles)
    
    set(handles.well_popup,'string',handles.obj.wellLabels)

else
    return
end

function handles=load_bc_data(hObject, eventdata, handles)
% extract barcode columns from the fcs file based on the barcode key

set(handles.parent,'pointer','watch')
drawnow

try
    handles.obj=handles.obj.find_bc_cols_by_mass;
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
sample_size=100000;
handles.obj=handles.obj.load_bcs(sample_size);

handles.obj=handles.obj.compute_debarcoding('bcs');

handles.obj=handles.obj.normalize_by_pop('bcs');

% handles.obj=handles.obj.normalize_bcs('bcs');
% %calculates normbcs from bcs

handles.obj=handles.obj.compute_debarcoding('normbcs');
%calculates bcinds from normbcs

%% 20140904 -- main cofactor update 
% handles.obj=handles.obj.calculate_cofactors; 
% 
% handles.obj=handles.obj.recofactor;
% %calculates cofactored_bcs from bcs and cofactors
% 
% handles.obj=handles.obj.normalize_bcs('cofactored_bcs');
% %calculates normbcs from cofactored_bcs
% 
% handles.obj=handles.obj.compute_debarcoding;
% %calculates bcind from normbcs
%% end 

% compute mahalanobis distances
handles.obj=handles.obj.compute_mahal;

handles.obj=handles.obj.compute_well_abundances;

% plot well abundances
set(handles.plottype,'SelectedObject',handles.separation_plot)
handles=separation_plot(handles);

drawnow
set(handles.parent,'pointer','arrow')
drawnow

set(handles.x_popup,'value',1)
set(handles.y_popup,'value',2)

if any(cellfun(@isempty,handles.obj.m(handles.obj.bc_cols)))
    set(handles.x_popup,'string',handles.obj.c(handles.obj.bc_cols))
    set(handles.y_popup,'string',handles.obj.c(handles.obj.bc_cols))
else
    set(handles.x_popup,'string',handles.obj.m(handles.obj.bc_cols))
    set(handles.y_popup,'string',handles.obj.m(handles.obj.bc_cols))
end

handles.leg=cell(1,handles.obj.num_masses);
for i=1:handles.obj.num_masses
    m_i=handles.obj.m{handles.obj.bc_cols(i)};
    if ~isempty(m_i)
        handles.leg{i}=m_i;
    else
        handles.leg{i}=handles.obj.c{handles.obj.bc_cols(i)};
    end
end


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

% update the value of the mahal cutoff stored in the scd object whenever it
% is changed in the GUI
handles.obj.mahal_cutoff_val=str2double(get(handles.mahal_cutoff,'string'));

within_cutoffs=handles.obj.mahal<handles.obj.mahal_cutoff_val & handles.obj.deltas>handles.obj.sep_cutoff;
for i=1:handles.obj.num_codes
    handles.obj.well_yield(i)=handles.obj.sample_ratio*nnz(handles.obj.bcind==i & within_cutoffs);
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

bar(yield_axis,1:handles.obj.num_codes,handles.obj.well_yield,'facecolor',[0 0.5 0.4])
ylabel('Cell count')
set(yield_axis,'xlim',[0 handles.obj.num_codes+1],'xtick',[])
yl=get(yield_axis,'ylim');
text(1:handles.obj.num_codes,-yl(2)/15*ones(1,handles.obj.num_codes),handles.obj.wellLabels,'rotation',315)


total_yield=round(100*sum(handles.obj.well_yield)/size(handles.obj.x,1));
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
if popval < length(handles.obj.wellLabels);
    set(handles.well_popup,'value',popval+1);
end
plot_button_Callback(hObject, eventdata, handles)


