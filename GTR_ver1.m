function varargout = GTR_ver1(varargin)
% GTR_ver1, Ground Thermal Respond program
% Created by Husni Firmansyah Sutrisno
% This program analyze the thermal respond of borehole with output of
% ground thermal conductivity (lambda) and borehole thermal resistance(Rb) 
%
% The input is excel file or other dataset in form of .txt file type. First
% page is for calculating the Thermal Response Test (TRT) of the borehole.
% Input parameter such as borehole diameter and active depth, heat
% injection starting & stopping time, and which time period should the 
% calculation take place. 
%
% Case selection:
% 
% Case 1: Optimize the lambda and Rb value at the same time
% Case 2: Calculate the Rb value only. Input estimating value of lambda
% Case 3: Calculate the lambda value only. Input value of Rb
% Case 4: See the temperature respond by selecting specific lambda and Rb
%         value.
% 
% DTRT Page
% Input file: Depth section, Temperature measured each section (TTMAT), and
% time data (time_hours). File type should be in .txt file.
% 
% Calculate DTS will initiate the calculation for each depth section.
% Control the result of each section on the table Lambda & Rb per depth.
% Curve fitting between temperature measured and responds will be shown on
% axes beside the table. The curve will update for each section result.
% 
% Result of Lambda and Rb per each depth section will be shown on graphs on
% result page tab. 
% 
% GTR_VER1 MATLAB code for GTR_ver1.fig
%       GTR_VER1, by itself, creates a new GTR_VER1 or raises the existing
%       singleton*.
% 
%      H = GTR_VER1 returns the handle to a new GTR_VER1 or the handle to
%      the existing singleton*.
%
%      GTR_VER1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GTR_VER1.M with the given input arguments.
%
%      GTR_VER1('Property','Value',...) creates a new GTR_VER1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GTR_ver1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GTR_ver1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 20-Sep-2018 10:52:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GTR_ver1_OpeningFcn, ...
                   'gui_OutputFcn',  @GTR_ver1_OutputFcn, ...
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

function GTR_ver1_OpeningFcn(hObject, eventdata, handles, varargin)

clc

% input variable (constant)

set(handles.timeoffset,'String','0')
set(handles.lambda_set,'Enable','off')
set(handles.rb_set,'Enable','off')
set(handles.firstrow,'Enable','off')
set(handles.lastrow,'Enable','off')
handles.offset = 0;
handles.case = 1;
handles.row_del = 0;
handles.output = hObject;

% initialise Logo/Icon
axes(handles.axes7)
imshow('iconGEO.PNG');
axes(handles.axes8)
imshow('BD_logo.PNG')

% Initialise tabs -- Do Not Delete
handles.tabManager = TabManager( hObject );
guidata(hObject, handles);

function varargout = GTR_ver1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function data_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to data_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
% table = get(handles.data_table,'Data');
% refreshDisplay(table,handles,1);


function plot_type_Callback(hObject, eventdata, handles)
% --- Set the plot based on the selection ---
% after browse TRT data, the plot is automatically showing temp profile
% Input is only the selection of the plot type

val = get(hObject,'Value');
strlist = get(hObject,'String');

tin = getappdata(handles.browse,'tin');
tout = getappdata(handles.browse,'tout');
flow = getappdata(handles.browse,'flow')*1000;
time = getappdata(handles.browse,'time');
power = getappdata(handles.browse,'power');
tmeas = getappdata(handles.browse,'tmeas');
cla reset
switch strlist{val};
    case 'Temp profile'        
        axes(handles.axes1)
        plot(time,tmeas,'.')
        title('Temp uppmätt')
        ylabel('Temperatur [C]')
        xlabel('Tid [tim]')
        legend('Temp uppmätt','Location','Southeast')
        tableData = [time,tmeas];
        set(handles.data_table,'data',tableData)
        set(handles.data_table,'ColumnName',{'Time';'Temp'})
            
    case 'Power & Flow profile'
        f_new = figure;
        f_new.Visible = 'off';
        plot(time,power,'.')
        ylim([min(power) max(power)*1.5])
        ylabel('Effekt [W]')
        yyaxis right
        plot(time,flow,'-r')
        ylabel('Flöde [l/s]')
        xlabel('Tid [tim]')
        set(gca,'YColor','k')
        legend('Effekt (beräk.)','Flöde','Location','Southeast')
        file = fullfile(handles.PathTRT,'Power');
        print(f_new,file,'-dmeta')
        delete(f_new)
        axes(handles.axes1)
        plot(time,power,'.')
        ylim([min(power) max(power)*1.5])
        ylabel('Effekt [W]')
        yyaxis right
        plot(time,flow,'-r')
        title('Effekt & Flöde profil')        
        legend('Effekt (beräk.)','Flöde','Location','Southeast')
        ylabel('Flöde [l/s]')
        xlabel('Tid [tim]')
        set(gca,'YColor','k')
        tableData = [time,power,flow];
        set(handles.data_table,'data',tableData)
        set(handles.data_table,'ColumnName',{'Time';'Power';'Flow'})
    
    case 'Tin & Tout profile' 
        f_new = figure;
        f_new.Visible = 'off';
        plot(time,tin,'.',time,tout,'.')
        ylabel('Köldbärartemperatur [°C]')
        xlabel('Tid [tim]')
        yyaxis right
        plot(time,handles.tamb,'.k')
        legend('Tin','Tut','Tluft','Location','Southeast')
        ylabel('Uteluftstemperatur [°C]')
        set(gca,'YColor','k') 
        file = fullfile(handles.PathTRT,'Temperaturer');
        print(f_new,file,'-dmeta')
        delete(f_new)
        axes(handles.axes1)
        plot(time,tin,'.',time,tout,'.')
        title('Temp in & Temp out')
        ylabel('Köldbärartemperatur [°C]')
        xlabel('Tid [tim]')
        yyaxis right
        plot(time,handles.tamb,'.k')
        legend('Tin','Tut','Tluft','Location','Southeast')
        ylabel('Uteluftstemperatur [°C]')
        set(gca,'YColor','k')
        tableData = [time,tin,tout];
        set(handles.data_table,'data',tableData)
        set(handles.data_table,'ColumnName',{'Time';'Tin';'Tout'})
     
    case 'Respons (Curve Fitting)'  
        try 
         DTM = getappdata(handles.calculate,'DTM');
         dT = getappdata(handles.calculate,'dT');
         axes(handles.axes1)
         plot(time,dT,'.',time,DTM,'.');
         ylabel('Temperaturskillnad mellan KB och ostörd mark [K]')
         xlabel('Tid [tim]')
         legend('DTuppmätt', 'Respons','Location','Southeast')
         title('Curve Fitting TRT')
         tableData = [time,DTM,dT];
         set(handles.data_table,'data',tableData)
         set(handles.data_table,'ColumnName',{'Time';'Response';'DT measured'})
        catch 
            warndlg('Please run the TRT analysis first','Warning guys here','Modal')
        end
end

handles.plot_type_val = val;
guidata(hObject,handles)

function diameter_Callback(hObject, eventdata, handles)
% --- Input the borehole diameter in meter unit ---
handles.r1 = str2double(get(hObject,'String'))/2;
guidata(hObject,handles);

function baseline_temp_Callback(hObject, eventdata, handles)
% --- Input the T0 temperature in celcius ---
T0 = str2double(get(hObject,'String'));
setappdata(hObject,'T0',T0)
set(handles.baseline_temp_dts,'String',T0)
guidata(hObject,handles)

function baseline_temp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function heat_start_Callback(hObject, eventdata, handles)
% --- Input the heat injection start in hour ---
heatstart = str2double(get(hObject,'String'));
setappdata(hObject,'heatstart',heatstart)
set(handles.heat_start_dts,'String',heatstart)
guidata(hObject,handles)

function heat_stop_Callback(hObject, eventdata, handles)
% --- Input the heat injection stop in hour ---
heatstop = str2double(get(hObject,'String'));
setappdata(hObject,'heatstop',heatstop)
set(handles.heat_stop_dts,'String',heatstop)
guidata(hObject,handles)

function opt_start_Callback(hObject, eventdata, handles)
% --- Input the start of optimization period for calculation in hour ---
optstart = str2double(get(hObject,'String'));
setappdata(hObject,'optstart',optstart)
set(handles.opt_start_dts,'String',optstart)
guidata(hObject,handles)

function opt_stop_Callback(hObject, eventdata, handles)
% --- Input the stop of optimization period for calculation in hour ---
optstop = str2double(get(hObject,'String'));
setappdata(hObject,'optstop',optstop)
set(handles.opt_stop_dts,'String',optstop)
guidata(hObject,handles)

function edit7_Callback(hObject, eventdata, handles)
% --- This box only for displaying the iteration of the TRT analysis ---
a = num2str(getappdata(handles.calculate,'lambda'));
b = num2str(getappdata(handles.calculate,'Rb'));
data = [a,b];
set(hObject,'String',data)

function depth_Callback(hObject, eventdata, handles)
% --- Input the depth of borehole in meter ---
handles.H = str2double(get(hObject,'String'));
guidata(hObject,handles)

function fluid_Callback(hObject, eventdata, handles)
% --- This function will calculate the power based on the coolprop function ---
% selecting ETHANOL WATER or ETHYLENE GLYCOL WATER will require input on
% the concentration box. Otherwise, WATER doesnot need such input, Leave
% blank on the Concentration.
% Once proper selection is chosen, the table of Data Set will update value
% with power each time. Notice the table change after selecting the fluid
% type properly

val_fluid = get(hObject,'Value');
strlist = get(hObject,'String');

tin = getappdata(handles.browse,'tin');
tout = getappdata(handles.browse,'tout');
flow = getappdata(handles.browse,'flow');
time = getappdata(handles.browse,'time');

tableData = get(handles.data_table,'data');

switch strlist{val_fluid};
    case 'Water'
        fluid='Water';
        power=power_fluid(tout,tin,flow,fluid);
        
    case 'Ethanol Water'
        try
            waitfor(handles.concentration,'String')
            a = get(handles.concentration,'String');
            fluid=['INCOMP::' 'MEA2-' a '%'];
            power=power_fluid(tout,tin,flow,fluid);
        
        catch ME
            disp(ME.message)
        end
        
    case 'Ethylene Glycol Water'
        try
            waitfor(handles.concentration,'String')
            a = get(handles.concentration,'String');
            fluid=['INCOMP::' 'MEG2-' a '%'];
            power=power_fluid(tout,tin,flow,fluid);
        
        catch ME
            disp(ME.message)
        end      
            
end
        if isempty(power) == 1
            power = getappdata(handles.browse,'power');
        end
        tableData = [time,power];
        set(handles.data_table,'data',tableData)
        set(handles.data_table,'ColumnName',{'Time';'Power'})
        setappdata(hObject,'power',power)
        assignin('base','power2',power)

function concentration_Callback(hObject, eventdata, handles)
% --- Input the concentration value of the fluid type used in percent ---
% for example, 25.66% Ethanol Water, type 25.66 in the box
concent = get(hObject,'String');
setappdata(hObject,'val',concent)
guidata(hObject,handles)

function calculate_Callback(hObject, eventdata, handles)
% --- This button start the TRT analysis calculation ---
% make sure every input is done and correct
% Things to check:
%               BROWSE DATA
%               CASE TYPE selection type
%               PARAMETER INPUTS
%               FLUID TYPE SELECTION
%               CONCENTRATION value as needed
% 
% If CASE TYPE selection is not LAMBDA & RB OPT, check the input value on 
% LAMBDA or RB box
% Delete the value inside LAMBDA and RB after each CASE TYPE or before
% running the calculation again.
% 
% The result of Lambda and Rb value is displayed on the LAMBDA and RB box
% on the right side of the page. The curve fitting between response and
% measured temperature is showed on the axes. Click the SAVE FIGURE to save
% the figure to the folder

tic
time = getappdata(handles.browse,'time');
Noptstart = getappdata(handles.opt_start,'optstart');
Noptstop = getappdata(handles.opt_stop,'optstop');
tmeas = getappdata(handles.browse,'tmeas');
power2 = getappdata(handles.fluid,'power');

if isempty(power2) == 1
    power2 = getappdata(handles.browse,'power');
end

T0 = getappdata(handles.baseline_temp,'T0');
heatstart = getappdata(handles.heat_start,'heatstart');
heatstop = getappdata(handles.heat_stop,'heatstop');

H = handles.H;
r1 = handles.r1;
%
heatstart=find(time<=heatstart,1,'last');
heatstop=find(time<=heatstop,1,'last');
Noptstart=find(time<=Noptstart,1,'last');
Noptstop=find(time<=Noptstop,1,'last');
%

    switch handles.case
        case 1
            f=msgbox('Optimize Lambda and Rb');
            [lambda, Rb, DTM, dT] = NR_v4(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1);
            delete(f)
            
        case 2
            f=msgbox('Optimize Rb only');
            lambda = handles.lambda_set_val;
            [Rb, DTM, dT] = case2(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1,lambda);
            lambda_end = Rb;
            lambda_end(:,:) = lambda;
            lambda = lambda_end;
            delete(f)
            
        case 3
            f=msgbox('Optimize Lambda only');
            Rb = handles.rb_set_val;
            if handles.recovery == 0
                [lambda, DTM, dT] = case3(time,Noptstart,Noptstop,heatstart,heatstop,tmeas,power2,T0,H,r1,Rb);
            else
                [lambda, DTM, dT] = case3_RP(time,Noptstart,Noptstop,heatstart,heatstop,tmeas,power2,T0,H,r1,Rb);
            end
%             [lambda, DTM, dT] = case3(time,Noptstart,Noptstop,heatstart,heatstop,tmeas,power2,T0,H,r1,Rb);
            Rb_end = lambda;
            Rb_end(:,:) = Rb;
            Rb = Rb_end;
            delete(f)
            
        case 4
            f=msgbox('Calculating respond only');
            lambda = handles.lambda_set_val;
            Rb = handles.rb_set_val;
            [DTM, dT] = case4(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1,lambda,Rb);
            delete(f)
        
    end

    
% --- RMSD
RMSD = sqrt(mean((DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop)).^2));
% RMSD = sqrt((sum((DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop)).^2)/(Noptstop-Noptstart-1)));
% RMSD = sqrt((DTM-dT).^2
set(handles.RMSD,'String',num2str(RMSD))
setappdata(hObject,'lambda',lambda);
setappdata(hObject,'Rb',Rb);
setappdata(hObject,'DTM',DTM);
setappdata(hObject,'dT',dT);
Table_iter = [lambda',Rb'];
set(handles.edit7,'data',Table_iter)
axes(handles.axes1)
plot(time,dT,'.',time,DTM,'.');
legend('DTuppmätt', 'Respons','Location','Southeast')

toc
msgbox(({'Calculation Completed';['Total Calculation Time = ' num2str(toc) ' seconds']}),'Success')
setappdata(hObject,'toc',toc)

set(handles.lambda_val,'String',num2str(lambda(end)));
set(handles.rb_val,'String',num2str(Rb(end)));

function browse_Callback(hObject, eventdata, handles)
% --- Browse the input data ---
% Valid input data type is excel in term of xlsx, xlsm
% please check if the total row data in excel is not more than 5000
% otherwise the code below have to be changed

msgbox('Please select the excel file type','modal')
pause(1)
    
[FileName,PathName] = uigetfile({'*.xlsm;*.csv;*.xlsx'},'TRT Input Data');
% if FileName == 0
%     msgbox('Please select the correct file','modal')
%     pause(2)
%     close
% return
% end
excel = fullfile(PathName,FileName);
% answer = questdlg('Do you need correction for the input data? (into each 10 minutes measurement)');
% if answer == 'Yes'
%     M = xlsread(excel);
%     [ii,jj]=size(M);
%     tid = 0:1/6:M(ii,1)+1/6;
%     N = NaN(length(tid),jj);
%     N(:,1) = tid';
%     for k = 2:jj
%         a1 = 1;
%         for i = 2:length(tid)
%             a2 = find(M(:,1)<=tid(i),1,'last');
%             N(i,k) = mean(M(a1:a2,k));
%             a1 = a2+1;
%         end 
%     end
%        
% end

handles.PathTRT = PathName;
msgbox('Loading data, please wait','modal');
data = readall(datastore(excel));
time = data.Time;
power = data.Power;
tmeas = data.Temp;
tin = data.Tin;
tout = data.Tout;
flow = data.Flow;
% Ft = xlsread(excel,'Data','A21:K10000'); % --- check the max row in the input file, otherwise change K5000 row
% time = Ft(:,1)-Ft(1); 
% power = Ft(:,2);
% tmeas = Ft(:,3);
% tin = Ft(:,4);
% tout = Ft(:,5);
% tamb = Ft(:,6);
% handles.tamb = tamb;
% flow = Ft(:,7)/1000;
close
msgbox('Loading data completed');
pause(0.5)
close

% --- Plot temp profile automatically
axes(handles.axes1)
plot(time,tmeas,'.')
title('Temp measured')
ylabel('Temp [C]')
xlabel('Time [hour]')
legend('Temp meas')
tableData = [time,tmeas];
set(handles.data_table,'data',tableData)
set(handles.data_table,'ColumnName',{'Time';'Temp'})

setappdata(hObject,'power',power)
setappdata(hObject,'time',time)
setappdata(hObject,'tmeas',tmeas)
setappdata(hObject,'tin',tin)
setappdata(hObject,'tout',tout)
setappdata(hObject,'flow',flow)
handles.tout = tout;
guidata(hObject,handles)

function browse_dtrt_Callback(hObject, eventdata, handles)
% --- Browse data input for DTRT analysis
% Valid file type is txt format
% First file input is Depth, then TTMAT and last timehours
% check the number of row of each input
% Depth row must be the same with TTMAT row
% Timehours row must be the same with TTMAT column

[file,path] = uigetfile('*.txt','Depth');
msgbox('Select Depth file in txt format','modal');
depth = load(fullfile(path,file));
uiwait(msgbox(['Total row (depth): ' num2str(length(depth))],'modal'))
[file,path] = uigetfile('*.txt','TTMAT');
TTMAT = load(fullfile(path,file));
uiwait(msgbox(['Total row (TTMAT): ' num2str(size(TTMAT))],'modal'))
[file,path] = uigetfile('*.txt','Time_hours');
time_hours = load(fullfile(path,file));
uiwait(msgbox(['Total row (time_hours): ' num2str(length(time_hours))],'modal'))

% --- Change the depth matrix from row to column only
[m,~]=size(depth);
if m<2
    depth = depth';
end

% --- Change the DTS timehours matrix from row to column only
[m,~]=size(time_hours);
if m<2
    time_hours = time_hours';
end

timehours = getappdata(handles.browse,'time');
tempcmedel = TTMAT;
tmeas = getappdata(handles.browse,'tmeas');

axes(handles.axes6)
surf(time_hours,depth,TTMAT,'EdgeColor','none')
colormap jet
colorbar
xlabel('Tid [tim]')
ylabel('Djup [m]')
zlabel('Temp [°C]')   
zlim([0 20])

axes(handles.dtrt_axes)
plot(timehours,tmeas,time_hours,tempcmedel(10,:))
legend('TRT','DTS')

set(handles.lambda_table,'data',depth);
set(handles.lambda_table,'ColumnName',{'Depth [m]'})
set(handles.text17,'String',['Number of section: ' num2str(length(depth))]);
set(handles.active_depth,'String',['Active depth: ' num2str(-depth(end)) ' m']);

setappdata(hObject,'timehours',timehours)
setappdata(hObject,'tempcmedel',tempcmedel)
setappdata(hObject,'depth',depth)
setappdata(hObject,'TTMAT',TTMAT)
setappdata(hObject,'time_hours',time_hours)

handles.tempcmedel = tempcmedel;
handles.depth = depth;
handles.timehours_dts = time_hours;
handles.timehours_trt = timehours;
guidata(hObject,handles)

function calculate_dts_Callback(hObject, eventdata, handles)
% --- Calculate the DTS based on DTRT input data --- %
%  Check the input value especially the number of row between Depth, TTMAT,
%  and time_hours from BROWSE DTRT INPUT
timehours_trt = getappdata(handles.browse,'time');
heatstart = handles.heatstart_dts;
heatstop = handles.heatstop_dts;
Noptstart = handles.optstart_dts;
Noptstop = handles.optstop_dts;
power2 = getappdata(handles.fluid,'power');
power2 = getappdata(handles.browse,'power');


T0 = getappdata(handles.t0_prof,'T0');
H = handles.H;
r1 = handles.r1;
depth = getappdata(handles.browse_dtrt,'depth');
tempcmedel = getappdata(handles.browse_dtrt,'tempcmedel');
time_hours_dtrt = getappdata(handles.browse_dtrt,'time_hours');
time = time_hours_dtrt;


heatstart=find(time<=heatstart,1,'last');
heatstop=find(time<=heatstop,1,'last');
Noptstart=find(time<=Noptstart,1,'last');
Noptstop=find(time<=Noptstop,1,'last');



[~,fpath]=uiputfile('*.bmp','Save Graph for each depth');

% --- Power2 interpolation
power2 = interp1(timehours_trt,power2,time_hours_dtrt);

% --- Specific section analysis
firstrow = 1;
lastrow = length(depth);
if handles.checkbox2.Value > 0
    firstrow = handles.firstrow_val;
    lastrow = handles.lastrow_val;
    set(handles.text17,'String',['Number of section: '...
        num2str(length(depth))])
    set(handles.active_depth,'String',['Active depth: '...
        num2str(depth(lastrow)-depth(firstrow)) ' m'])
end

% --- calculation time begin
tic
lambda_end = NaN(length(depth),1);
Rb_end=lambda_end;
Nrow = length(tempcmedel(1,:));
DTM=zeros(Nrow,1);
% tableData = zeros(length(depth),3);

% ---  waitbar
f = waitbar(0,'0','Name','Calculating DTS',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

for k=firstrow:lastrow
        
tmeas = tempcmedel(k,:).';

if length(T0) < 2
    T0 = getappdata(handles.baseline_temp_dts,'T0');
else
    T0 = T0(k,1);
end 

disp(['Depth  = ', num2str(depth(k))])



%
    switch handles.case
        case 1
            
            [lambda, Rb, DTM, dT] = NR_v4(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1);
            
            
        case 2
           
            lambda = handles.lambda_set_val;
            [Rb, DTM, dT] = case2(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1,lambda);
            lambda_end = Rb;
            lambda_end(:,:) = lambda;
            lambda = lambda_end;
            
            
        case 3
            
            Rb = handles.rb_set_val;
            power2 = getappdata(handles.browse,'power');
            [lambda, DTM, dT] = case3(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1,Rb);
                        
            lambda_end(k,1)=lambda(end);
            Rb_end(k,1)=Rb(end);
            
        case 4
            
            lambda = handles.lambda_set_val;
            Rb = handles.rb_set_val;
            [DTM, dT] = case4(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1,lambda,Rb);
            
        
    end
%








% [lambda, Rb, DTM, dT] = NR_v4(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1);

lambda_end(k,1)=lambda(end);
Rb_end(k,1)=Rb(end);

T0 = getappdata(handles.t0_prof,'T0'); % --- Reset the T0 value

    if getappdata(f,'canceling')
        f1 = findall(0,'type','figure','tag','TMWWaitbar');
        delete(f1)
        break        
    end
    waitbar(k/length(depth),f,{'Depth '; num2str(depth(k))})


% --- set table update
tableData(k,:) = [depth(k,1),lambda_end(k,1),Rb_end(k,1)];
set(handles.lambda_table,'data',tableData)
set(handles.lambda_table,'ColumnName',{'Depth';'Lambda';'Rb'})

%%% set axes update
axes(handles.dtrt_axes)
plot(time_hours_dtrt,DTM,time_hours_dtrt+handles.offset,dT);
ylim([round(min(dT),1,'significant')-1 round(max(dT),1,'significant')+1])
legend('Response', 'DT measured','Location','southeast')
title(['Depth = ' num2str(depth(k)) ' Lambda = ' num2str(lambda_end(k))...
    ' Rb = ' num2str(Rb_end(k))])

% --- create GIF
% drawnow
% frame = getframe(handles.dtrt_axes);
% im{k} = frame2im(frame);
% filename = 'All_section.gif'; % Specify the output file name
% [A,map] = rgb2ind(im{k},256);
%     if k == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.5);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);
%     end

% --- Plot and save each depth/section
fname = ['Depth',num2str(depth(k)),'.bmp'];
file = fullfile(fpath,fname);
ax_old = gca;
f_new = figure;
f_new.Visible = 'off';
ax_new = copyobj(ax_old,f_new);
set(ax_new,'Position','default')
title(['Depth' num2str(depth(k)) ' lambda = ' num2str(lambda_end(k))...
    ' Rb = ' num2str(Rb_end(k))])
legend({'Response', 'DT measured'},'Location','southeast')
print(f_new,file,'-dpng')
delete(f_new)

% --- Plot Lambda and Rb for every depth
axes(handles.lambda_axes)
plot(lambda_end,depth,'.')
xlabel('Lambda')
ylabel('Depth (m)')
xlim([1 6])
ylim([min(depth) 0])
axes(handles.rb_axes)
plot(Rb_end,depth,'.')
xlabel('Rb')
ylabel('Depth (m)')
xlim([0 0.2])
ylim([min(depth) 0])

end
f1 = findall(0,'type','figure','tag','TMWWaitbar'); % -- delete the waitbar
delete(f1)
toc

msgbox(({'Calculation Completed';['Total Calculation Time = ' num2str(toc/60)...
    ' minutes']}),'Success')

lambda_end(isnan(lambda_end)) = [];
set(handles.lambda_ave,'String',num2str(mean(lambda_end)))

Rb_end(isnan(Rb_end)) = [];
set(handles.rb_ave,'String',num2str(mean(Rb_end)))

setappdata(hObject,'lambda_dts',lambda_end);
setappdata(hObject,'Rb_dts',Rb_end);
setappdata(hObject,'DTM_dts',DTM);
setappdata(hObject,'dT_dts',dT);

assignin('base','lambda',lambda_end)
assignin('base','Rb',Rb_end)
assignin('base','depth',depth)
guidata(hObject,handles)

function baseline_temp_dts_Callback(hObject, eventdata, handles)
% --- Input the T0 temperature in celcius or select the temp profile ---
T0 = str2double(get(hObject,'String'));
setappdata(hObject,'T0',T0)
set(handles.baseline_temp_dts,'String',T0)
guidata(hObject,handles)

function heat_start_dts_Callback(hObject, eventdata, handles)
% --- Input the heat injection start or select on the graph ---
handles.heatstart_dts = str2double(get(hObject,'String'));
guidata(hObject,handles)

function heat_stop_dts_Callback(hObject, eventdata, handles)
% --- Input the heat injection stop or select on the graph ---
handles.heatstop_dts = str2double(get(hObject,'String'));
guidata(hObject,handles)

function opt_start_dts_Callback(hObject, eventdata, handles)
% --- Input the start time of optimization period or select on the graph ---
handles.optstart_dts = str2double(get(hObject,'String'));
guidata(hObject,handles)

function opt_stop_dts_Callback(hObject, eventdata, handles)
% --- Input the stop time of optimization period or select on the graph ---
handles.optstop_dts = str2double(get(hObject,'String'));
guidata(hObject,handles)

function timeoffset_Callback(hObject, eventdata, handles)
% --- Input the time offset between TRT and DTRT ---
% check on the lower graph if any offset between both. input the value in
% hour and the graph will change based on the value input.
handles.offset = str2double(get(hObject,'String'));
timehours = getappdata(handles.browse,'time');
time_hours = getappdata(handles.browse_dtrt,'time_hours');
tmeas = getappdata(handles.browse,'tmeas');
tempcmedel = getappdata(handles.browse_dtrt,'tempcmedel');
axes(handles.dtrt_axes)
plot(timehours-handles.offset,tmeas,time_hours,tempcmedel(10,:))
xlim([0 max(time_hours)])
legend('TRT','DTS')
guidata(hObject,handles)

function lambda_set_Callback(hObject, eventdata, handles)
% --- Lambda value if the case selected for TRT analysis is RB OPT or 
% RESPOND ONLY 
handles.lambda_set_val = str2num(get(hObject,'String'));
guidata(hObject,handles)

function rb_set_Callback(hObject, eventdata, handles)
% --- Rb value if the case selected for TRT analysis is LAMBDA OPT or
% RESPOND ONLY
handles.rb_set_val = str2num(get(hObject,'String'));
guidata(hObject,handles)

function case1_Callback(hObject, eventdata, handles)
% --- Case for LAMBDA & RB OPT ---
set(handles.lambda_set,'Enable','off')
set(handles.rb_set,'Enable','off')
handles.case = 1;
guidata(hObject,handles)

function case2_Callback(hObject, eventdata, handles)
% --- Case for RB OPT ---
set(handles.rb_set,'Enable','off')
set(handles.lambda_set,'Enable','on')
handles.case = 2;
guidata(hObject,handles)

function case3_Callback(hObject, eventdata, handles)
% --- Case for LAMBDA OPT ---
set(handles.lambda_set,'Enable','off')
set(handles.rb_set,'Enable','on')
handles.case = 3;
guidata(hObject,handles)

function case4_Callback(hObject, eventdata, handles)
% --- Case for RESPOND ONLY ---
set(handles.lambda_set,'Enable','on')
set(handles.rb_set,'Enable','on')
handles.case = 4;
guidata(hObject,handles)

function savelambda_Callback(hObject, eventdata, handles)
% --- Button to save the lambda graph for DTRT
[fname,fpath]=uiputfile('*.emf');
file = fullfile(fpath,fname);
axes(handles.lambda_axes)
ax_old = gca;
f_new = figure;
f_new.Visible = 'off';
ax_new = copyobj(ax_old,f_new);
% title('Lambda per depth')
set(ax_new,'OuterPosition','default','XAxisLocation','top')
grid on
print(f_new,file,'-dmeta')
delete(f_new)


function saverb_Callback(hObject, eventdata, handles)
% --- Button to save the lambda graph for DTRT
[fname,fpath]=uiputfile('*.emf');
file = fullfile(fpath,fname);
axes(handles.rb_axes)
ax_old = gca;
f_new = figure;
f_new.Visible = 'off';
ax_new = copyobj(ax_old,f_new);
set(ax_new,'OuterPosition','default','XAxisLocation','top')
grid on
% title('Rb per depth')
print(f_new,file,'-dmeta')
delete(f_new)

function savetrt_Callback(hObject, eventdata, handles)
% --- Button to save the curve fitting graph for TRT analysis
time = getappdata(handles.browse,'time');
tmeas = getappdata(handles.browse,'tmeas');
DTM = getappdata(handles.calculate,'DTM');
dT = getappdata(handles.calculate,'dT');
[fname,fpath]=uiputfile('*.emf');
file = fullfile(fpath,fname);
f_new = figure('OuterPosition',[100 100 1200 600]);
f_new.Visible = 'off';

switch handles.plot_type_val
    case 2
        plot(time,tmeas,'.')
        ylabel('Temperatur [C]')
        xlabel('Tid [tim]')
        legend('Temp uppmätt','Location','Southeast')      
        print(f_new,file,'-dmeta')
        uiwait(msgbox({'Saved on' file},'modal'))
    case 3
        uiwait(msgbox({'Saved on' handles.PathTRT 'Power.png'},'modal'))
    case 4
        uiwait(msgbox({'Saved on' handles.PathTRT 'Temperaturer.png'},'modal'))
    case 5
        plot(time,dT,'.',time,DTM,'.');
        ylabel('Temperaturskillnad mellan KB och ostörd mark [K]')
        xlabel('Tid [tim]')
        legend('DTuppmätt', 'Respons','Location','Best')  
        print(f_new,file,'-dmeta')
        uiwait(msgbox({'Saved on' file},'modal'))
end

delete(f_new)

function [lambda, Rb, DTM, dT] = NR_v4(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1)
% ---------
% This code is specifically functioning for Matlab version of 2016 or lower 
% due to matrix operation is not compatible with older version than 2017
% for newer version please choose NR_v3
% ---------

gamma = 0.5772; 
rho_cp = 2.16*10^6;
alpha = 1.72*10^-6; 
time = time*3600;
dT = tmeas-T0;
dP=(power2(2:end)-power2(1:end-1))/H;  
Nrow = size(time,1);
DTM=zeros(Nrow,1);


N=Noptstop-Noptstart+1;                

time2=time-time(heatstart); 
logtime2=log(time2);
power2_mean=mean(power2(Noptstart:Noptstop));

reg = [ones(N,1) logtime2(Noptstart:Noptstop)];  
bm = reg\dT(Noptstart:Noptstop);
b = bm(1);
m = bm(2);
lambda(1) = power2_mean/(H*4*pi*m);
a = lambda(1)/rho_cp;                       
Rb(1) = H*(b-m*(log(4*a/r1^2)-gamma))/power2_mean;

if lambda(1) <= 1 || lambda(1) == Inf 
        lambda(1) = 3;
        Rb(1) = 0.08;
end

if Rb(1) <= 0.005
    Rb(1) = 0.08;
end

Ltol = 0.0001;
Rtol = 0.0001;
count = 1;
df2dr = sum(2*power2(Noptstart:Noptstop).^2/(H^2));

timeMj=bsxfun(@times,ones(Noptstop-1,N),time(1:Noptstop-1)); 
timeMi=bsxfun(@times,ones(Noptstop-1,N),time(Noptstart:Noptstop)'); 
timeM=timeMi-timeMj;
timeM(timeM<0)=0;

while count==1 || abs(lambda(count)-lambda(count-1))/lambda(count) > Ltol || abs(Rb(count)-Rb(count-1))/Rb(count) > Rtol
    
    alpha = lambda(count)/rho_cp;
    disp(['Lambda = ', num2str(lambda(count)),'    Rb = ', num2str(Rb(count))])
        
    DTMM=bsxfun(@times,dP(1:Noptstop-1),(1/(4*pi*lambda(count)).*expint(r1^2./(4*alpha*timeM))+Rb(count)));
    
    DTMM(isnan(DTMM))=0;
    DTM(Noptstart:Noptstop)=sum(DTMM)';
    x=r1^2./(4*alpha*timeM); 
    d=expint(x);
    d(isnan(d))=0;
    e=exp(-x);
    e(isnan(e))=0;
    c=bsxfun(@times,dP(1:Noptstop-1),(e-d));
    c=sum(c);
    g=bsxfun(@times,dP(1:Noptstop-1),(bsxfun(@times,(x-3),e)+2*d));
    g(isnan(g))=0;
    g=sum(g);
    f1=c'.*(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop));
    f1=sum(f1);
    f2=power2(Noptstart:Noptstop).*(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop))/H;
    f2=sum(f2);
    df1dl=bsxfun(@times,(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop)),g')+c'.^2/(4*pi*lambda(count));
    df1dl=sum(df1dl);    
    df1dr = bsxfun(@times,c',power2(Noptstart:Noptstop)/H);
    df1dr=sum(df1dr);
    f1 = f1/(2*pi*lambda(count)^2);
    f2 = 2*f2;
    df1dl = df1dl / (2*pi*lambda(count)^3);
    df1dr = df1dr/(2*pi*lambda(count)^2);
    det=(df1dl*df2dr-df1dr^2);
    lambda(count+1)=min(max(lambda(count)-(f1*df2dr-f2*df1dr)/det,1),10);
    Rb(count+1)=min(max(Rb(count)-(f2*df1dl-f1*df1dr)/det,0),0.3);
%     [a A]=min(lambda(count+1)-lambda);
%     [b B]=min(Rb(count+1)-Rb);
%     if A==B && abs(a)/lambda(A)<Ltol && abs(b)/lambda(B)<Rtol && count>4 && A~=count
%         lambda(count+1)=lambda(1)*(0.5+rand(1));
%         Rb(count+1)=Rb(1)*(0.5+rand(1));
%     end
    count=count+1;  
    if count-1> 19
        if strcmp(questdlg(['the algorithm does not seem to be converging for this section, do you want' ...
                ' to skip it?'],'YesNo'),'Yes')
            lambda(count)=lambda(count-1);
            Rb(count)=Rb(count-1);
        else
            lambda(count+1)=lambda(1)*(0.5+rand(1));
            Rb(count+1)=Rb(1)*(0.5+rand(1));
        end
    end
end

timeMj=bsxfun(@times,ones(Noptstart-2,Noptstart-1),time(1:Noptstart-2));
timeMi=bsxfun(@times,ones(Noptstart-2,Noptstart-1),time(1:Noptstart-1)');
timeM=timeMi-timeMj;
timeM(timeM<0)=0;
DTMM=bsxfun(@times,dP(1:Noptstart-2),(1/(4*pi*lambda(count)).*expint(r1^2./(4*alpha*timeM))+Rb(count)));
DTMM(isnan(DTMM))=0;
DTM(1:Noptstart-1)=sum(DTMM)';

timeMj=bsxfun(@times,ones(Nrow-1,Nrow-Noptstop+1),time(1:Nrow-1));
timeMi=bsxfun(@times,ones(Nrow-1,Nrow-Noptstop+1),time(Noptstop:Nrow)');
timeM=timeMi-timeMj;
timeM(timeM<0)=0;
DTMM=bsxfun(@times,dP(1:Nrow-1),(1/(4*pi*lambda(count)).*expint(r1^2./(4*alpha*timeM))+Rb(count)));
DTMM(isnan(DTMM))=0;
DTM(Noptstop:Nrow)=sum(DTMM)';

function [Rb, DTM, dT] = case2(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1,lambda)
% --------- Calculating Rb only -------------------------------------------
% This code is specifically functioning for Matlab version of 2016 or lower 
% due matrix operation is not compatible with older version than 2017
% for newer version please choose NR_v3
% ---------

gamma = 0.5772; 
rho_cp = 2.16*10^6;
alpha = 1.72*10^-6; 
time = time*3600;
dT = tmeas-T0;
dP=(power2(2:end)-power2(1:end-1))/H;  
Nrow = size(time,1);
DTM=zeros(Nrow,1);

N=Noptstop-Noptstart+1;                 

time2=time-time(heatstart); 
logtime2=log(time2);
power2_mean=mean(power2(Noptstart:Noptstop));

reg = [ones(N,1) logtime2(Noptstart:Noptstop)]; 
bm = reg\dT(Noptstart:Noptstop);
b = bm(1);
m = bm(2);
a = lambda/rho_cp;                     
Rb(1) = H*(b-m*(log(4*a/r1^2)-gamma))/power2_mean;

if Rb(1) <= 0
        
        Rb(1) = 0.08;
end

Ltol = 0.001;
Rtol = 0.001;
count = 1;
df2dr = sum(2*power2(Noptstart:Noptstop).^2/(H^2));

timeMj=bsxfun(@times,ones(Noptstop-1,N),time(1:Noptstop-1)); 
timeMi=bsxfun(@times,ones(Noptstop-1,N),time(Noptstart:Noptstop)'); 
timeM=timeMi-timeMj;
timeM(timeM<0)=0;

while count==1 || abs(Rb(count)-Rb(count-1))/Rb(count) > Rtol
    
    alpha = lambda/rho_cp;
    disp(['Lambda = ', num2str(lambda),'    Rb = ', num2str(Rb(count))]) 
    
    DTMM=bsxfun(@times,dP(1:Noptstop-1),(1/(4*pi*lambda).*expint(r1^2./(4*alpha*timeM))+Rb(count)));
    
    DTMM(isnan(DTMM))=0;
    DTM(Noptstart:Noptstop)=sum(DTMM)';
    x=r1^2./(4*alpha*timeM); 
    d=expint(x);
    d(isnan(d))=0;
    e=exp(-x);
    e(isnan(e))=0;
    c=bsxfun(@times,dP(1:Noptstop-1),(e-d));
    c=sum(c);
    g=bsxfun(@times,dP(1:Noptstop-1),(bsxfun(@times,(x-3),e)+2*d));
    g(isnan(g))=0;
    g=sum(g);
    f1=c'.*(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop));
    f1=sum(f1);
    f2=power2(Noptstart:Noptstop).*(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop))/H;
    f2=sum(f2);
    df1dl=bsxfun(@times,(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop)),g')+c'.^2/(4*pi*lambda);
    df1dl=sum(df1dl);    
    df1dr = bsxfun(@times,c',power2(Noptstart:Noptstop)/H);
    df1dr=sum(df1dr);
    f1 = f1/(2*pi*lambda^2);
    f2 = 2*f2;
    df1dl = df1dl / (2*pi*lambda^3);
    df1dr = df1dr/(2*pi*lambda^2);
    det=(df1dl*df2dr-df1dr^2);
    Rb(count+1)=min(max(Rb(count)-(f2*df1dl-f1*df1dr)/det,0.005),0.3);

    count=count+1;  
    if count-1> 19
        if strcmp(questdlg(['the algorithm does not seem to be converging for this section, do you want' ...
                ' to skip it?'],'YesNo'),'Yes')
            Rb(count)=Rb(count-1);
        else
            Rb(count+1)=Rb(1)*(0.5+rand(1));
        end
    end
end

timeMj=bsxfun(@times,ones(Noptstart-2,Noptstart-1),time(1:Noptstart-2));
timeMi=bsxfun(@times,ones(Noptstart-2,Noptstart-1),time(1:Noptstart-1)');
timeM=timeMi-timeMj;
timeM(timeM<0)=0;
DTMM=bsxfun(@times,dP(1:Noptstart-2),(1/(4*pi*lambda).*expint(r1^2./(4*alpha*timeM))+Rb(count)));
DTMM(isnan(DTMM))=0;
DTM(1:Noptstart-1)=sum(DTMM)';

timeMj=bsxfun(@times,ones(Nrow-1,Nrow-Noptstop+1),time(1:Nrow-1));
timeMi=bsxfun(@times,ones(Nrow-1,Nrow-Noptstop+1),time(Noptstop:Nrow)');
timeM=timeMi-timeMj;
timeM(timeM<0)=0;
DTMM=bsxfun(@times,dP(1:Nrow-1),(1/(4*pi*lambda).*expint(r1^2./(4*alpha*timeM))+Rb(count)));
DTMM(isnan(DTMM))=0;
DTM(Noptstop:Nrow)=sum(DTMM)';

    
function [lambda, DTM, dT] = case3(time,Noptstart,Noptstop,heatstart,heatstop,tmeas,power2,T0,H,r1,Rb)
% --------- Calculating Lambda only ---------------------------------------
% This code is specifically functioning for Matlab version of 2016 or lower 
% due matrix operation is not compatible with older version than 2017
% for newer version please choose NR_v3
% ---------

rho_cp = 2.16*10^6;
alpha = 1.72*10^-6; 
time = time*3600;
dT = tmeas-T0;
dP =(power2(2:end)-power2(1:end-1))/H;  
Nrow = size(time,1);
DTM=zeros(Nrow,1);

N=Noptstop-Noptstart+1;                 

time2=time-time(heatstart); 
logtime2=log(time2);
power2_mean=mean(power2(Noptstart:Noptstop));

reg = [ones(N,1) logtime2(Noptstart:Noptstop)];  
bm = reg\dT(Noptstart:Noptstop);
b = bm(1);
m = bm(2);
lambda(1) = power2_mean/(H*4*pi*m);
          
if lambda(1) <= 0 || lambda(1) == Inf
        lambda(1) = 3;  
end

Ltol = 0.0001;

count = 1;
df2dr = sum(2*power2(Noptstart:Noptstop).^2/(H^2));

timeMj=bsxfun(@times,ones(Noptstop-1,N),time(1:Noptstop-1)); 
timeMi=bsxfun(@times,ones(Noptstop-1,N),time(Noptstart:Noptstop)'); 
timeM=timeMi-timeMj;
timeM(timeM<0)=0;

while count==1 || abs(lambda(count)-lambda(count-1))/lambda(count) > Ltol 
    
    alpha = lambda(count)/rho_cp;
    disp(['Lambda = ', num2str(lambda(count)),'    Rb = ', num2str(Rb)])
  
    DTMM=bsxfun(@times,dP(1:Noptstop-1),(1/(4*pi*lambda(count)).*expint(r1^2./(4*alpha*timeM))+Rb));
    
    DTMM(isnan(DTMM))=0;
    DTM(Noptstart:Noptstop)=sum(DTMM)';
    x=r1^2./(4*alpha*timeM); %
    d=expint(x);
    d(isnan(d))=0;
    e=exp(-x);
    e(isnan(e))=0;
    c=bsxfun(@times,dP(1:Noptstop-1),(e-d));
    c=sum(c);
    g=bsxfun(@times,dP(1:Noptstop-1),(bsxfun(@times,(x-3),e)+2*d));
    g(isnan(g))=0;
    g=sum(g);
    f1=c'.*(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop));
    f1=sum(f1);
    f2=power2(Noptstart:Noptstop).*(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop))/H;
    f2=sum(f2);
    df1dl=bsxfun(@times,(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop)),g')+c'.^2/(4*pi*lambda(count));
    df1dl=sum(df1dl);    
    df1dr = bsxfun(@times,c',power2(Noptstart:Noptstop)/H);
    df1dr=sum(df1dr);
    f1 = f1/(2*pi*lambda(count)^2);
    f2 = 2*f2;
    df1dl = df1dl / (2*pi*lambda(count)^3);
    df1dr = df1dr/(2*pi*lambda(count)^2);
    det=(df1dl*df2dr-df1dr^2);
    lambda(count+1)=min(max(lambda(count)-(f1*df2dr-f2*df1dr)/det,1),10);
    
    count=count+1;  
    if count-1> 19
        if strcmp(questdlg(['the algorithm does not seem to be converging for this section, do you want' ...
                ' to skip it?'],'YesNo'),'Yes')
            lambda(count)=lambda(count-1);
        else
            lambda(count+1)=lambda(1)*(0.5+rand(1));
        end
    end
end

timeMj=bsxfun(@times,ones(Noptstart-2,Noptstart-1),time(1:Noptstart-2));
timeMi=bsxfun(@times,ones(Noptstart-2,Noptstart-1),time(1:Noptstart-1)');
timeM=timeMi-timeMj;
timeM(timeM<0)=0;
DTMM=bsxfun(@times,dP(1:Noptstart-2),(1/(4*pi*lambda(count)).*expint(r1^2./(4*alpha*timeM))+Rb));
DTMM(isnan(DTMM))=0;
DTM(1:Noptstart-1)=sum(DTMM)';

timeMj=bsxfun(@times,ones(Nrow-1,Nrow-Noptstop+1),time(1:Nrow-1));
timeMi=bsxfun(@times,ones(Nrow-1,Nrow-Noptstop+1),time(Noptstop:Nrow)');
timeM=timeMi-timeMj;
timeM(timeM<0)=0;
DTMM=bsxfun(@times,dP(1:Nrow-1),(1/(4*pi*lambda(count)).*expint(r1^2./(4*alpha*timeM))+Rb));
DTMM(isnan(DTMM))=0;
DTM(Noptstop:Nrow)=sum(DTMM)';

function [lambda, DTM, dT] = case3_RP(time,Noptstart,Noptstop,heatstart,heatstop,tmeas,power2,T0,H,r1,Rb)
% --------- Calculating Lambda only ---------------------------------------
% This code is specifically functioning for Matlab version of 2016 or lower 
% due matrix operation is not compatible with older version than 2017
% for newer version please choose NR_v3
% ---------

rho_cp = 2.16*10^6;
alpha = 1.72*10^-6; 
time = time*3600;
dT = tmeas-T0;
dP =(power2(2:end)-power2(1:end-1))/H;  
Nrow = size(time,1);
DTM=zeros(Nrow,1);

N=Noptstop-Noptstart+1;                 

time2=time-time(heatstart); 
logtime2=log(time2);
power2_mean=mean(power2(heatstart:heatstop));

reg = [ones(N,1) logtime2(Noptstart:Noptstop)];  
bm = reg\dT(Noptstart:Noptstop);
b = bm(1);
m = bm(2);
lambda(1) = power2_mean/(H*4*pi*m);
          
if lambda(1) <= 0 || lambda(1) == Inf
        lambda(1) = 3;  
end

Ltol = 0.0001;

count = 1;
df2dr = sum(2*power2_mean.^2/(H^2));

timeMj=bsxfun(@times,ones(Noptstop-1,N),time(1:Noptstop-1)); 
timeMi=bsxfun(@times,ones(Noptstop-1,N),time(Noptstart:Noptstop)'); 
timeM=timeMi-timeMj;
timeM(timeM<0)=0;

while count==1 || abs(lambda(count)-lambda(count-1))/lambda(count) > Ltol 
    
    alpha = lambda(count)/rho_cp;
    disp(['Lambda = ', num2str(lambda(count)),'    Rb = ', num2str(Rb)])
  
    DTMM=bsxfun(@times,dP(1:Noptstop-1),(1/(4*pi*lambda(count)).*expint(r1^2./(4*alpha*timeM))+Rb));
    
    DTMM(isnan(DTMM))=0;
    DTM(Noptstart:Noptstop)=sum(DTMM)';
    x=r1^2./(4*alpha*timeM); %
    d=expint(x);
    d(isnan(d))=0;
    e=exp(-x);
    e(isnan(e))=0;
    c=bsxfun(@times,dP(1:Noptstop-1),(e-d));
    c=sum(c);
    g=bsxfun(@times,dP(1:Noptstop-1),(bsxfun(@times,(x-3),e)+2*d));
    g(isnan(g))=0;
    g=sum(g);
    f1=c'.*(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop));
    f1=sum(f1);
    f2=power2_mean.*(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop))/H;
    f2=sum(f2);
    df1dl=bsxfun(@times,(DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop)),g')+c'.^2/(4*pi*lambda(count));
    df1dl=sum(df1dl);    
    df1dr = bsxfun(@times,c',power2_mean/H);
    df1dr=sum(df1dr);
    f1 = f1/(2*pi*lambda(count)^2);
    f2 = 2*f2;
    df1dl = df1dl / (2*pi*lambda(count)^3);
    df1dr = df1dr/(2*pi*lambda(count)^2);
    det=(df1dl*df2dr-df1dr^2);
    lambda(count+1)=min(max(lambda(count)-(f1*df2dr-f2*df1dr)/det,1),10);
    
    count=count+1;  
    if count-1> 19
        if strcmp(questdlg(['the algorithm does not seem to be converging for this section, do you want' ...
                ' to skip it?'],'YesNo'),'Yes')
            lambda(count)=lambda(count-1);
        else
            lambda(count+1)=lambda(1)*(0.5+rand(1));
        end
    end
end

timeMj=bsxfun(@times,ones(Noptstart-2,Noptstart-1),time(1:Noptstart-2));
timeMi=bsxfun(@times,ones(Noptstart-2,Noptstart-1),time(1:Noptstart-1)');
timeM=timeMi-timeMj;
timeM(timeM<0)=0;
DTMM=bsxfun(@times,dP(1:Noptstart-2),(1/(4*pi*lambda(count)).*expint(r1^2./(4*alpha*timeM))+Rb));
DTMM(isnan(DTMM))=0;
DTM(1:Noptstart-1)=sum(DTMM)';

timeMj=bsxfun(@times,ones(Nrow-1,Nrow-Noptstop+1),time(1:Nrow-1));
timeMi=bsxfun(@times,ones(Nrow-1,Nrow-Noptstop+1),time(Noptstop:Nrow)');
timeM=timeMi-timeMj;
timeM(timeM<0)=0;
DTMM=bsxfun(@times,dP(1:Nrow-1),(1/(4*pi*lambda(count)).*expint(r1^2./(4*alpha*timeM))+Rb));
DTMM(isnan(DTMM))=0;
DTM(Noptstop:Nrow)=sum(DTMM)';

function [DTM, dT] = case4(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1,lambda,Rb)
% --------- Calculating Respond Only provided by Lambda and Rb set --------
% This code is specifically functioning for Matlab version of 2016 or lower 
% due matrix operation is not compatible with older version than 2017
% for newer version please choose NR_v3
% ---------

rho_cp = 2.16*10^6;
time = time*3600;
dT = tmeas-T0;
dP=(power2(2:end)-power2(1:end-1))/H;  
Nrow = size(time,1);
DTM=zeros(Nrow,1);

N=Noptstop-Noptstart+1;                 

timeMj=bsxfun(@times,ones(Noptstop-1,N),time(1:Noptstop-1)); 
timeMi=bsxfun(@times,ones(Noptstop-1,N),time(Noptstart:Noptstop)'); 
timeM=timeMi-timeMj;
timeM(timeM<0)=0;

alpha = lambda/rho_cp;
disp(['Lambda = ', num2str(lambda),'    Rb = ', num2str(Rb)])

DTMM=bsxfun(@times,dP(1:Noptstop-1),(1/(4*pi*lambda).*expint(r1^2./(4*alpha*timeM))+Rb));
DTMM(isnan(DTMM))=0;
DTM(Noptstart:Noptstop)=sum(DTMM)';

timeMj=bsxfun(@times,ones(Noptstart-2,Noptstart-1),time(1:Noptstart-2));
timeMi=bsxfun(@times,ones(Noptstart-2,Noptstart-1),time(1:Noptstart-1)');
timeM=timeMi-timeMj;
timeM(timeM<0)=0;
DTMM=bsxfun(@times,dP(1:Noptstart-2),(1/(4*pi*lambda).*expint(r1^2./(4*alpha*timeM))+Rb));
DTMM(isnan(DTMM))=0;
DTM(1:Noptstart-1)=sum(DTMM)';

timeMj=bsxfun(@times,ones(Nrow-1,Nrow-Noptstop+1),time(1:Nrow-1));
timeMi=bsxfun(@times,ones(Nrow-1,Nrow-Noptstop+1),time(Noptstop:Nrow)');
timeM=timeMi-timeMj;
timeM(timeM<0)=0;
DTMM=bsxfun(@times,dP(1:Nrow-1),(1/(4*pi*lambda).*expint(r1^2./(4*alpha*timeM))+Rb));
DTMM(isnan(DTMM))=0;
DTM(Noptstop:Nrow)=sum(DTMM)';



function power2=power_fluid(tout,tin,flow,fluid)
% --- Calculating Power using CoolProp property ---------------------------
% Change the path "CoolPpath" on the new computer

    CoolPpath='C:\Users\Carry\Documents\Bengt_Dahlgren\GUI\Coolprop_matlab';

    
    a=1;
while a==1
    try
        addpath(CoolPpath);
        CoolProp.PropsSI('D', 'T', 293.15, 'P', 101325,'Water');
        a=0;
    catch ME
        waitfor(msgbox(['You may have entered a wrong path for CoolProp.' char(10) ...
        ['Check that the Coolprop path is ' CoolPpath] char(10) ...
        'If you have not installed CoolProp, visit:' char(10) 'http://www.coolprop.org/coolprop/wrappers/MATLAB/index.html']));
        if strcmp(questdlg('Do you want to enter a new path for CoolProp?' ...
        ,'Temperature profile updated','YesNo'),'Yes')
            CoolPpath=uigetdir('C:\Program Files\CoolProp_Matlab','Select the CoolProp directory');
        else
            return
        end
    end
end
    delta_T = (tout-tin);
    
    power2 = CoolProp.PropsSI('D', 'T', mean(mean([tin tout]))+273.15, 'P', 101325,fluid)*...
        CoolProp.PropsSI('C', 'T', mean(mean([tin tout]))+273.15, 'P', 101325,fluid).*delta_T.*flow;

    power2(1,1)=0;
    
function radiobutton5_Callback(hObject, eventdata, handles)
% --- Selecting the heat start value from graph
if get(hObject,'Value')==1
    [x,~] = ginput(1);
    set(handles.heat_start_dts,'String',num2str(x))
    axes(handles.dtrt_axes)
    ax = gca;
    y = ax.YLim;
    line([x x],y,'Color','c');    
    handles.heatstart_dts = x;
end
 guidata(hObject,handles)

function radiobutton6_Callback(hObject, eventdata, handles)
% --- Selecting the heat stop value from graph
if get(hObject,'Value')==1
    [x,~] = ginput(1);
    set(handles.heat_stop_dts,'String',num2str(x))
    axes(handles.dtrt_axes)
    ax = gca;
    y = ax.YLim;
    line([x x],y,'Color','y'); 
    handles.heatstop_dts = x;
end
 guidata(hObject,handles)
 
function radiobutton7_Callback(hObject, eventdata, handles)
% --- Selecting the opt start value from graph
if get(hObject,'Value')==1
    [x,~] = ginput(1);
    set(handles.opt_start_dts,'String',num2str(x))
    axes(handles.dtrt_axes)
    ax = gca;
    y = ax.YLim;
    line([x x],y,'Color','k');
    handles.optstart_dts = x;
end
 guidata(hObject,handles)


function radiobutton8_Callback(hObject, eventdata, handles)
% --- Selecting the opt stop value from graph
if get(hObject,'Value')==1
    [x,~] = ginput(1);
    set(handles.opt_stop_dts,'String',num2str(x))
    axes(handles.dtrt_axes)
    ax = gca;
    y = ax.YLim;
    line([x x],y,'Color','k');   
    handles.optstop_dts = x;
end
 guidata(hObject,handles)

function checkbox1_Callback(hObject, eventdata, handles)
% --- button to delete all line produced when selecting the heat start,
% heat stop, opt start, opt stop, value by clicking on the graph directly
F = findall(0,'type','Line','DisplayName','');
delete(F)
f1 = findall(0,'type','figure','tag','TMWWaitbar'); % -- delete the waitbar
delete(f1)

function checkbox2_Callback(hObject, eventdata, handles)
% --- Checkbox to allow access to box of delete row
if get(hObject,'Value')==1
    set(handles.firstrow,'Enable','on')
    set(handles.lastrow,'Enable','on')
    handles.selectrow = 1;
else
    set(handles.firstrow,'Enable','off')
    set(handles.lastrow,'Enable','off')
end
guidata(hObject,handles)

function firstrow_Callback(hObject, eventdata, handles)
handles.firstrow_val = str2double(get(hObject,'String'));
guidata(hObject,handles)

function lastrow_Callback(hObject, eventdata, handles)
handles.lastrow_val = str2double(get(hObject,'String'));
depth = getappdata(handles.browse_dtrt,'depth');
TTMAT = getappdata(handles.browse_dtrt,'TTMAT');
time_hours = getappdata(handles.browse_dtrt,'time_hours');

depth = depth(handles.firstrow_val:handles.lastrow_val);
TTMAT = TTMAT(handles.firstrow_val:handles.lastrow_val,:);
set(handles.lambda_table,'data',depth);
set(handles.text17,'String',['Number of section: ' num2str(length(depth))])
set(handles.active_depth,'String',['Active depth: ' num2str(abs(depth(end)-depth(1))) ' m'])
axes(handles.axes6)
surf(time_hours,depth,TTMAT,'EdgeColor','none')
colormap jet
colorbar
xlabel('Tid [tim]')
ylabel('Djup [m]')
zlabel('Temp [ºC]')
% zlim([0 40])         
guidata(hObject,handles)

function row_delete_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in t0_prof.
function t0_prof_Callback(hObject, eventdata, handles)
% --- selecting the temp profile for baseline temp (T0) ---
%  click on the graph to look the profile on that specific time.
if get(hObject,'Value')==1
    plot_trt_dts_Callback(handles.plot_trt_dts, eventdata, handles)
    [x,~] = ginput(1);
    disp(['Temp profile on  ' num2str(x) '  hours'])
    depth = getappdata(handles.browse_dtrt,'depth');
    time_hours = getappdata(handles.browse_dtrt,'time_hours');
    tempcmedel = getappdata(handles.browse_dtrt,'tempcmedel');
    a = max(find(time_hours<x));
    axes(handles.dtrt_axes)
    plot(tempcmedel(:,a),depth);
    T0 = tempcmedel(:,a);
    setappdata(hObject,'T0',T0);
    set(handles.baseline_temp_dts,'String','Profile selected')
end

% --- Executes on button press in plot_trt_dts.
function plot_trt_dts_Callback(hObject, eventdata, handles)
% --- Click button to get the graph of DTS & TRT on the DTRT page 
tmeas = getappdata(handles.browse,'tmeas');
timehours = getappdata(handles.browse,'time');
time_hours = getappdata(handles.browse_dtrt,'time_hours');
tempcmedel = getappdata(handles.browse_dtrt,'tempcmedel');
axes(handles.dtrt_axes)
plot(timehours-handles.offset,handles.tout,'.',time_hours,tempcmedel(handles.firstrow_val,:),'.')
legend('TRT','DTS')

function quit_Callback(hObject, eventdata, handles)
% --- Quit button on TRT Page
close(ancestor(hObject,'figure'))

function pushbutton9_Callback(hObject, eventdata, handles)
% --- Quit button on DTRT Page
close(ancestor(hObject,'figure'))


function timeoffset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_stop_dts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_start_dts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function heat_stop_dts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function heat_start_dts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function baseline_temp_dts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lambda_val_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rb_val_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function concentration_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fluid_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function depth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function opt_stop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.optstop_dts = str2double(get(hObject,'String'));
guidata(hObject,handles)

function opt_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.optstart_dts = str2double(get(hObject,'String'));
guidata(hObject,handles)

function heat_stop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.heatstop_dts = str2double(get(hObject,'String'));
guidata(hObject,handles)

function heat_start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.heatstart_dts = str2double(get(hObject,'String'));
guidata(hObject,handles)

function diameter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plot_type_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rb_set_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lambda_set_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function data_table_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function calculate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function firstrow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lastrow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lambda_ave_Callback(hObject, eventdata, handles)




function lambda_ave_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rb_ave_Callback(hObject, eventdata, handles)


function rb_ave_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RMSD_Callback(hObject, eventdata, handles)

% Hints: get(hObject,'String') returns contents of RMSD as text
%        str2double(get(hObject,'String')) returns contents of RMSD as a double


% --- Executes during object creation, after setting all properties.
function RMSD_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start1_Callback(hObject, eventdata, handles)
handles.Noptstart1 = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function start1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start2_Callback(hObject, eventdata, handles)
handles.Noptstart2 = str2double(get(hObject,'String'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function start2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function stop21_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function stop11_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function stop11_Callback(hObject, eventdata, handles)
handles.Noptstop1 = str2double(get(hObject,'String'));
guidata(hObject,handles)

function stop12_Callback(hObject, eventdata, handles)
handles.Noptstop12 = str2double(get(hObject,'String'));
guidata(hObject,handles)

function stop21_Callback(hObject, eventdata, handles)
handles.Noptstop2 = str2double(get(hObject,'String'));
guidata(hObject,handles)

function stop22_Callback(hObject, eventdata, handles)
handles.Noptstop22 = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes on button press in calculate_sensi.
function calculate_sensi_Callback(hObject, eventdata, handles)
Noptstop1 = handles.Noptstop1:5:handles.Noptstop12; 
Noptstop2 = handles.Noptstop2:5:handles.Noptstop22;

period1 = Noptstop1-handles.Noptstart1;
period2 = Noptstop2-handles.Noptstart2;

time = getappdata(handles.browse,'time');
tmeas = getappdata(handles.browse,'tmeas');
power2 = getappdata(handles.fluid,'power');
T0 = getappdata(handles.baseline_temp,'T0');
heatstart = getappdata(handles.heat_start,'heatstart');
H = handles.H;
r1 = handles.r1;

lambda1 = ones(length(Noptstop1),1);
Rb1 = lambda1;
lambda2 = ones(length(Noptstop2),1);
Rb2 = lambda2;

for i = 1:length(Noptstop1)
    Noptstart = handles.Noptstart1;
    Noptstop = Noptstop1(i);
    [lambda, Rb] = NR_v4(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1);
    lambda1(i) = lambda(end);
    Rb1(i) = Rb(end);
end

% for i = 1:length(Noptstop2)
%     Noptstart = handles.Noptstart2;
%     Noptstop = Noptstop2(i);
%     [lambda, Rb] = NR_v4(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1);
%     lambda2(i) = lambda(end);
%     Rb2(i) = Rb(end);
% end
x = {['Starttid ' num2str(handles.Noptstart1)]};
y = {num2str(handles.Noptstart2)};
limit_x = [(period1(1)-10) (period1(end)+10)];

% --- save figure automatically
f_new = figure;
f_new.Visible = 'off';
subplot(2,1,1)
plot(period1,lambda1,'+') %,period2,lambda2,'x')
ylabel('\lambda [W/mK]')
xlabel('Optimeringstid [tim]')
% legend(x,'Location','best')     %num2str(handles.Noptstart1),num2str(handles.Noptstart2))
% xlim([0 max(period1)+20])
xlim(limit_x)
ylim([1 6])
grid on
% file = fullfile(handles.PathTRT,'Sensitivity Lambda');
% print(f_new,file,'-dmeta')
% delete(f_new)

% f_new = figure;
% f_new.Visible = 'off';
subplot(2,1,2)
plot(period1,Rb1,'+')   %,period2,Rb2,'x')
ylabel('Rb [Km/W]')
xlabel('Optimeringstid [tim]')
% legend(x)                 %num2str(handles.Noptstart1),num2str(handles.Noptstart2))
% xlim([0 max(period1)+20])
xlim(limit_x)
ylim([0 0.3])
grid on
file = fullfile(handles.PathTRT,'Sensitivity');
print(f_new,file,'-dmeta')
delete(f_new)
% ---


axes(handles.axes11)
plot(period1,lambda1,'+')               %,period2,lambda2,'x')
ylabel('\lambda [W/mK]')
xlabel('Optimeringstid [tim]')
% legend(x)                         %num2str(handles.Noptstart1),num2str(handles.Noptstart2))
% xlim([0 max(period1)+20])
xlim(limit_x)
ylim([1 6])
grid on
axes(handles.axes12)
plot(period1,Rb1,'+')                   %,period2,Rb2,'x')
ylabel('Rb [Km/W]')
xlabel('Optimeringstid [tim]')
% legend(x)                                 %num2str(handles.Noptstart1),num2str(handles.Noptstart2))
% xlim([0 max(period1)+20])
xlim(limit_x)
ylim([0 0.3])
grid on

% --- Executes during object creation, after setting all properties.
function stop12_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function stop22_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save3d.
function save3d_Callback(hObject, eventdata, handles)
[fname,fpath]=uiputfile('*.png');
file = fullfile(fpath,fname);
axes(handles.axes6)
ax_old = gca;
f_new = figure;
f_new.Visible = 'off';
ax_new = copyobj(ax_old,f_new);
c = colorbar;
colormap jet
% c.Label.String = 'Temp';
% c.Label.Rotation = 0;
% c.Limits = [0 15];
set(ax_new,'OuterPosition','default')
grid on
print(f_new,file,'-dpng')
delete(f_new)

function ostord_Callback(hObject, eventdata, handles)
handles.ostord = str2double(get(hObject,'String'));

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ostord_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function cirkulation_Callback(hObject, eventdata, handles)
handles.cirkulation = str2double(get(hObject,'String'));

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function cirkulation_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function varme1_Callback(hObject, eventdata, handles)
handles.varme1 = str2double(get(hObject,'String'));

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function varme1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function varme2_Callback(hObject, eventdata, handles)
handles.varme2 = str2double(get(hObject,'String'));

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function varme2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ater1_Callback(hObject, eventdata, handles)
handles.ater1 = str2double(get(hObject,'String'));

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ater1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ater2_Callback(hObject, eventdata, handles)
handles.ater2 = str2double(get(hObject,'String'));

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function ater2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in calculate_temp.
function calculate_temp_Callback(hObject, eventdata, handles)
depth = handles.depth;
temp = handles.tempcmedel;
starting_time = datenum(handles.date_start);
formatOut = 'yyyy-mm-dd HH:MM';
l1 = datestr(handles.ostord*60*60*1.1574e-05+starting_time,formatOut);
l2 = datestr(handles.cirkulation*60*60*1.1574e-05+starting_time,formatOut);
l3 = datestr(handles.varme1*60*60*1.1574e-05+starting_time,formatOut);
l4 = datestr(handles.varme2*60*60*1.1574e-05+starting_time,formatOut);
l5 = datestr(handles.ater1*60*60*1.1574e-05+starting_time,formatOut);
l6 = datestr(handles.ater2*60*60*1.1574e-05+starting_time,formatOut);


ostord=find(handles.timehours_dts<=handles.ostord,1,'last');
cirkulation=find(handles.timehours_dts<=handles.cirkulation,1,'last');
varme1=find(handles.timehours_dts<=handles.varme1,1,'last');
varme2=find(handles.timehours_dts<=handles.varme2,1,'last');
ater1=find(handles.timehours_dts<=handles.ater1,1,'last');
ater2=find(handles.timehours_dts<=handles.ater2,1,'last');

% --- save figure automatically
f_new = figure;
f_new.Visible = 'off';
cla reset
hold on
plot(temp(:,ostord),depth,'b')
plot(temp(:,cirkulation),depth,':b','LineWidth',1.5)
plot(temp(:,varme1),depth,'r')
plot(temp(:,varme2),depth,'-.r')
% plot(temp(:,ater1),depth,'k')
% plot(temp(:,ater2),depth,'--k')
% legend('ostörd','cirkulation','värme 1','värme 2','återhämtning 1','återhämtning 2','Location','Southeast')
legend({l1,l2,l3,l4},'Location','southeast')
xlabel('Temperatur [°C]')
ylabel('Djup [m]')
set(gca,'XAxisLocation','top');
% xlim([7 16])
hold off
file = fullfile(handles.PathTRT,'Temperature_profiles');
print(f_new,file,'-dmeta')
delete(f_new)
% ---

axes(handles.axes13);
cla reset
hold on
plot(temp(:,ostord),depth,'b')
plot(temp(:,cirkulation),depth,':b','LineWidth',1.5)
plot(temp(:,varme1),depth,'r')
plot(temp(:,varme2),depth,'-.r')
% plot(temp(:,ater1),depth,'k')
% plot(temp(:,ater2),depth,'--k')

% legend('ostörd','cirkulation','värme 1','värme 2','återhämtning 1','återhämtning 2','Location','Southeast')
legend(l1,l2,l3,l4,'Location','southeast')
xlabel('Temperatur [°C]')
ylabel('Djup [m]')
set(gca,'XAxisLocation','top');
% xlim([7 16])



function starting_time_Callback(hObject, eventdata, handles)
handles.date_start = get(hObject,'String');
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function starting_time_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in correct_dts.
function correct_dts_Callback(hObject, eventdata, handles)
row = handles.firstrow_val;
tempC_in = handles.tempcmedel(row,:)';
tout = interp1(handles.timehours_trt,handles.tout,handles.timehours_dts);
dev1 = tout - tempC_in;
[r, col] = size(handles.tempcmedel);
TTMAT_korr = zeros(size(handles.tempcmedel));
for i = 1:r
    TTMAT_korr(i,:)= handles.tempcmedel(i,:)+dev1(i);
end
dlmwrite('TTMAT_korr.txt',TTMAT_korr,'delimiter','\t')

function lambda_1_Callback(hObject, eventdata, handles)

function lambda_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lambda_2_Callback(hObject, eventdata, handles)

function lambda_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Rb_1_Callback(hObject, eventdata, handles)

function Rb_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Rb_2_Callback(hObject, eventdata, handles)

function Rb_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NN_Callback(hObject, eventdata, handles)

function NN_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function calculate_rmsd_Callback(hObject, eventdata, handles)
NN = str2double(get(handles.NN,'String'));
lambda_1 = str2double(get(handles.lambda_1,'String'));
lambda_2 = str2double(get(handles.lambda_2,'String'));
Rb_1 = str2double(get(handles.Rb_1,'String'));
Rb_2 = str2double(get(handles.Rb_2,'String'));
lambda_rmsd=linspace(lambda_1,lambda_2,NN);
Rb_rmsd=linspace(Rb_1,Rb_2,NN);
RMSD=zeros(NN);
time = getappdata(handles.browse,'time');
Noptstart = getappdata(handles.opt_start,'optstart');
Noptstop = getappdata(handles.opt_stop,'optstop');
tmeas = getappdata(handles.browse,'tmeas');
power2 = getappdata(handles.fluid,'power');
if isempty(power2) == 1
    power2 = getappdata(handles.browse,'power');
end

T0 = getappdata(handles.baseline_temp,'T0');
heatstart = getappdata(handles.heat_start,'heatstart');
H = handles.H;
r1 = handles.r1;
%
heatstart=find(time<=heatstart,1,'last');
Noptstart=find(time<=Noptstart,1,'last');
Noptstop=find(time<=Noptstop,1,'last');
%
tic
f=msgbox('Calculating RMSD factor');

for i = 1:NN;
    for j = 1:NN;
        lambda = lambda_rmsd(i);
        Rb = Rb_rmsd(j);
        [DTM, dT] = case4(time,Noptstart,Noptstop,heatstart,tmeas,power2,T0,H,r1,lambda,Rb);
        RMSD(i,j) = sqrt((sum((DTM(Noptstart:Noptstop)-dT(Noptstart:Noptstop)).^2)/(Noptstop-Noptstart-1)));
    end
end

axes(handles.axes_rmsd)
contourf(Rb_rmsd,lambda_rmsd,RMSD,100);
xlabel('Rb')
ylabel('\lambda')
colormap jet
colorbar
delete(f)
toc
msgbox(({'Calculation Completed';['Total Calculation Time = ' num2str(toc/60) ' minutes']}),'Finished')

% --- save figure automatically
f_new = figure;
f_new.Visible = 'off';
contourf(Rb_rmsd,lambda_rmsd,RMSD,100);
xlabel('Rb')
ylabel('\lambda')
colormap jet
colorbar
file = fullfile(handles.PathTRT,'RMSD');
print(f_new,file,'-dmeta')
delete(f_new)
% ---


% --- Executes on button press in recovery period selection.
function recovery_Callback(hObject, eventdata, handles)

if get(hObject,'Value')==1
    handles.recovery = 1;
else
    handles.recovery = 0;
end
guidata(hObject,handles)

