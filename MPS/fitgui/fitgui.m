function varargout = fitgui(varargin)
% Version 2 - 15.10.2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      Created by Ariel Nause for the use of lab A
%                    Version 2 created by Adiel Meyer 
%                          Tel Aviv University 2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changes in version 2:
%  - Error calculations are more accurate, by using the correlations between
%     the parameters.
%   - Covariance is computed and represented. 
%   - P-probability is computed and represented.
%   - Ndf is represented. (Number of Degrees of Freedom )
%   - Analytical error calculation of a straight line is calculated and
%     represented. (Errors are calculated by ignoring the errors in the X axis).
%   - Values and  errors are displayed in the same  line.
% 
%
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @fitgui_OpeningFcn, ...
    'gui_OutputFcn',  @fitgui_OutputFcn, ...
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



% --- Executes just before fitgui is made visible.
function fitgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% Choose default command line output for fitgui
handles.output = hObject;
handles.Data = FitGuiDataClass();

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = fitgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = handles.output;


function filename_edit_Callback(hObject, eventdata, handles)
filename_input = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function filename_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)


function load_edit_Callback(hObject, eventdata, handles)
% loads the data and saves the variables
%filename_input = get(handles.filename_edit,'String'); % gets the name of the data .dat/.txt file
[filename, pathname] = uigetfile( ...
{'*.txt;','Text Files (*.txt)';
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file');
 
filename_full = fullfile(pathname,filename);   
set(handles.filename_edit,'String', filename_full);
filename_input = get(handles.filename_edit,'String'); % gets the name of the data .dat/.txt file
data_given = load(filename_input);
if size(data_given,2)==4
    x = data_given(:,1); %set the range
    dx = data_given(:,2); % x errors
    y = data_given(:,3);  % Assigns the second column of mydata to a vector called 'y'
    dy = data_given(:,4); % y errors
end
if size(data_given,2)==2
    x = data_given(:,1);  % Assigns the first column of mydata to a vector called 'x'
    dx = zeros(length(data_given), 1); % x errors
    y = data_given(:,2);  % Assigns the second column of mydata to a vector called 'y'
    dy = zeros(length(data_given), 1); % y errors
end
if (size(data_given,2)~=2)&&(size(data_given,2)~=4)
    'Data file error'
end

handles.Data.x = x;
handles.Data.y = y;
handles.Data.dx = dx;
handles.Data.dy = dy;
guidata(hObject, handles)


function plotdata_pushbutton_Callback(hObject, eventdata, handles)
%Plots the data only - performs no fit at all. input as x dx y dy.

% clf(1);
if get(handles.PlotNewWindow,'Value');
    figure()
else
    figure(1);
end

data = handles.Data;
x = data.x; y = data.y; dx = data.dx; dy = data.dy;
%ploterr(x,y,dx,dy,'b.');
hold on
ploterr(x,y,dx,dy,'b.','hhxy',0.0001);
hold off
guidata(hObject, handles)


% --- Executes on button press in non_lin_checkbox.
function non_lin_checkbox_Callback(hObject, eventdata, handles)
%checkboxStatus = 0, if the box is unchecked,
%checkboxStatus = 1, if the box is checked
checkboxStatus = get(handles.non_lin_checkbox,'Value');
if(checkboxStatus)
    %if box is checked needs to get a function name
    set(handles.function_edit,'Visible' , 'on');
    set(handles.a0_edit,'Visible' , 'on');
    set(handles.text6,'Visible' , 'on');
	set(handles.text4,'Visible' , 'on');
    set(handles.PlotInitialPushButton,'Visible' , 'on');
    
else
    %if box is unchecked, no need for a function name
    set(handles.function_edit,'Visible', 'off');
    set(handles.a0_edit,'Visible' , 'off');
    set(handles.text6,'Visible' , 'off');
	set(handles.text4,'Visible' , 'off');
	set(handles.PlotInitialPushButton,'Visible' , 'off');
end
guidata(hObject, handles)

% --- Executes on button press in PlotInitialPushButton.
function PlotInitialPushButton_Callback(hObject, eventdata, handles)
% hObject    handle to PlotInitialPushButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function_input = get(handles.function_edit,'String'); % gets the name of the function .m file
a0 = str2num(get(handles.a0_edit,'String')); % this is the initial values guess as entered by user

data = handles.Data;
x = data.x; y = data.y; dx = data.dx; dy = data.dy;

checkboxStatusRange = get(handles.change_range_check_box,'Value');
if(checkboxStatusRange)
    %if box is checked needs to get a new x range
    x_min = str2num(get(handles.x_min_edit,'string')); % gets the x minimal value from user
    x_max = str2num(get(handles.x_max_edit,'string')); % gets the x maximal value from user
else
    %if box is unchecked, no need to get a new x range
    x_min = min(x); % gets the x minimal value from vector
    x_max = max(x); % gets the x maximal value from vector
end
ind = find( x >= x_min & x <= x_max ); 
x_fit = x(ind); %x_fit is the range the fit is performed upon

if min(x_fit) > 0
    plot_min = min(x_fit)*0.9;
else
    plot_min = min(x_fit)*1.1;
end
if max(x_fit) > 0
    plot_max = max(x_fit)*1.1;
else
    plot_max = max(x_fit)*0.9;
end
x_plot = linspace(plot_min , plot_max , 500);    

y_plot = feval(function_input,x_plot,a0);
%if get(handles.PlotNewWindow,'Value');
%    figure()
%else
%   figure(1);
%end
figure(100)
hold on
plot(x_plot,y_plot,'g-');        % Plots  the fitted curve
ploterr(x,y,dx,dy,'b.','hhxy',0.0001); % Plots also the data
hold off



function fit_pushbutton_Callback(hObject, eventdata, handles)
% performs the fitting and present the results - Not plotting!
filename_input = get(handles.filename_edit,'String'); % gets the name of the data .dat/.txt file
function_input = get(handles.function_edit,'String'); % gets the name of the function .m file
a0 = str2num(get(handles.a0_edit,'String')); % this is the initial values guess as entered by user
type = get(handles.non_lin_checkbox,'Value');

data = handles.Data;
x = data.x; y = data.y; dx = data.dx; dy = data.dy;

checkboxStatusRange = get(handles.change_range_check_box,'Value');
if(checkboxStatusRange)
    %if box is checked needs to get a new x range
    x_min = str2num(get(handles.x_min_edit,'string')); % gets the x minimal value from user
    x_max = str2num(get(handles.x_max_edit,'string')); % gets the x maximal value from user
else
    %if box is unchecked, no need to get a new x range
    x_min = min(x); % gets the x minimal value from vector
    x_max = max(x); % gets the x maximal value from vector
end
ind = find( x >= x_min & x <= x_max ); 
x_fit = x(ind); %x_fit is the range the fit is performed upon
dx_fit = dx(ind); % x_fit errors
y_fit = y(ind);  % the y values in the required range
dy_fit = dy(ind); % y fit errors

if min(x_fit) > 0
    plot_min = min(x_fit)*0.9;
else
    plot_min = min(x_fit)*1.1;
end
if max(x_fit) > 0
    plot_max = max(x_fit)*1.1;
else
    plot_max = max(x_fit)*0.9;
end
x_fit_plot = linspace(plot_min , plot_max , 500);    


if type == 0
	function_input ='lin';
    %x_fit_plot = x_fit; % for plot compatability with non linear fit
    [ainit,aerr,cov,chisq,y_fit_plot] = fitlin(x_fit,y_fit,dy_fit); % operates the chi^2 function fit
    aerr_closed=aerr;
    cov_closed=cov;
    a0=ainit;
	[a,aerr,cov,chisq,y_fit_plot] = fitnonlin(x_fit,x_fit_plot,y_fit,dx_fit,dy_fit,function_input,a0);
    RChiSquare = chisq/(length(x_fit)-length(a)) ;     %Reduced Chi-Square of your fit.
    % print results of fit to screen:
	
    fprintf('\n')
    disp(['initial parameters'' values:'])
    for x= 1:length(a0)
       fprintf('a%d = %f ± %f ', x, a0(x),aerr_closed(x));
    end
    fprintf('\n')
    disp(['fitted parameters'' values:'])
   for x= 1:length(a0)
       disp(['a' num2str(x) ' = ' , num2str(a(x) , 7),' ± ' num2str(aerr(x) , 7)])
   end
    if length(a0)> 1 
        for x= 1:length(a0)- 1
                for y =  x+1:length(a0)
                    disp(['cov(a' num2str(x) ',a', num2str(y), ') = ' , num2str(cov(x,y) , 7)])
                end
         end
    end
    disp(['chi^2 = ' , num2str(chisq , 5)])
    disp(['ndf = ' , num2str(length(x_fit)-length(a) , 5)])
    disp(['chi^2_reduced = ' , num2str(RChiSquare , 5)])
    disp(['p probability = ' , num2str(1 - chi2cdf(chisq,length(x_fit)-length(a)) , 5)])
end
if type == 1
    [a,aerr,cov,chisq,y_fit_plot] = fitnonlin(x_fit,x_fit_plot,y_fit,dx_fit,dy_fit,function_input,a0);    
    RChiSquare = chisq/(length(x_fit)-length(a)) ;     %Reduced Chi-Square of your fit.
    % print results of fit to screen:
    fprintf('\n')
    disp(['initial parameters'' values:'])
   for x= 1:length(a0)
       fprintf('a%d = %f  ', x, a0(x));
   end
   fprintf('\n')
    disp(['fitted parameters'' values:'])
   for x= 1:length(a0)
       disp(['a' num2str(x) ' = ' , num2str(a(x) , 7),' ± ' num2str(aerr(x) , 7)])
   end
       if length(a0)> 1 
        for x= 1:length(a0)- 1
                for y =  x+1:length(a0)
                    disp(['cov(a' num2str(x) ',a', num2str(y), ') = ' , num2str(cov(x,y) , 7)])
                end
         end
    end
    disp(['chi^2 = ' , num2str(chisq , 5)])
    disp(['ndf = ' , num2str(length(x_fit)-length(a) , 5)])
    disp(['chi^2_reduced = ' , num2str(RChiSquare , 5)])
    disp(['p probability = ' , num2str(1 - chi2cdf(chisq,length(x_fit)-length(a)) , 5)])
end

handles.Data.Fit.x_fit = x_fit;
handles.Data.Fit.y_fit = y_fit;
handles.Data.Fit.a = a;
handles.Data.Fit.aerr = aerr;
handles.Data.Fit.chisq = chisq;
handles.Data.Fit.RChiSquare = RChiSquare;
handles.Data.Fit.x_fit_plot = x_fit_plot;
handles.Data.Fit.y_fit_plot = y_fit_plot;

guidata(hObject, handles)

function [interpreterType] = GetInterpreterString(hObject, handles, popupName)
value = get(handles.(popupName), 'Value');
switch value
    case 2
        interpreterType = 'tex';
    case 3
        interpreterType = 'latex';
    otherwise
        interpreterType = 'none';
end

% --- Executes on button press in plot_pushbutton.
function plot_pushbutton_Callback(hObject, eventdata, handles)
% Plots the fit on the plotted data figure. Plots only the limited section
% x range the user enters

data = handles.Data;
x = data.x; y = data.y; dx = data.dx; dy = data.dy;

fit_data = data.Fit;

if get(handles.PlotNewWindow,'Value');
    figure()
else
    figure(1);
end
% hold on
% ploterr(x,y,dx,dy,'b.','hhxy',0.0001);
hold on
plot(fit_data.x_fit_plot,fit_data.y_fit_plot,'r-');        % Plots  the fitted curve
hold off
% axis([x(1),x(length(x)),min(y)-2*max(dy),max(y)+2*max(dy)]);
% axis([x_fit(1),x_fit(length(x_fit)),min(y_fit)-2*max(dy_fit),max(y_fit)+2*max(dy_fit)]);
% sets visible range of the plot

main_title = get(handles.main_title_edit,'String'); % gets the plot title
x_title = get(handles.x_title_edit,'String'); % gets the x title
y_title = get(handles.y_title_edit,'String'); % gets the y title

title(main_title,'fontsize',12, 'interpreter', GetInterpreterString(hObject, handles, 'titleInterpreterPopup')); % Places the title on the graph
xlabel(x_title ,'fontsize',12, 'interpreter', GetInterpreterString(hObject, handles, 'xTitleInterpreterPopup'));  % Labels the 'x' axis
ylabel(y_title ,'fontsize',12, 'interpreter', GetInterpreterString(hObject, handles, 'yTitleInterpreterPopup'));  % Labels the 'y' axis

% displays RChiSquare on graph
RChiSquare = fit_data.chisq/(length(fit_data.x_fit)-length(fit_data.a)) ;     %Reduced Chi-Square of your fit.
% legend( ['\chi^2_{RED} = ' , num2str(RChiSquare , 5)] , 'Location' , 'Best');
guidata(hObject, handles);


% --- Executes on button press in PlotResiduals.
function PlotResiduals_Callback(hObject, eventdata, handles)
% hObject    handle to PlotResiduals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = handles.Data;
x = data.x; y = data.y; dx = data.dx; dy = data.dy;

fit_data = data.Fit;

if get(handles.PlotNewWindow,'Value');
    figure()
else
    figure(200);
end

type = get(handles.non_lin_checkbox,'Value');
if type == 0  % linear fit
    y_res = y - feval('lin',x,fit_data.a);
elseif type == 1  % non-linear fit
    function_input = get(handles.function_edit,'String'); % gets the name of the function .m file
    y_res = y - feval(function_input,x,fit_data.a);
end

hold on
plot(fit_data.x_fit_plot,zeros(size(fit_data.y_fit_plot)),'r-');        % Plots a line of zeros
ploterr(x,y_res,dx,dy,'b.','hhxy',0.0001);  % Plots residuals
hold off

main_title = get(handles.main_title_edit,'String'); % gets the plot title
x_title = get(handles.x_title_edit,'String'); % gets the x title
y_title = get(handles.Residuals_y_title_edit,'String'); % gets the Residuals' y title

title([main_title , ' - Residuals Plot'] ,'fontsize',12, 'interpreter', GetInterpreterString(hObject, handles, 'titleInterpreterPopup')); % Places the title on the graph
xlabel(x_title ,'fontsize',12, 'interpreter', GetInterpreterString(hObject, handles, 'xTitleInterpreterPopup'));  % Labels the 'x' axis
ylabel(y_title ,'fontsize',12, 'interpreter', GetInterpreterString(hObject, handles, 'ResyTitleInterperterPopup'));  % Labels the 'y' axis
%%%%%%%%%%%%%%%%%%%%%  Set Axes & titles   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main_title_edit_Callback(hObject, eventdata, handles)
% This function creates the main title from user
guidata(hObject, handles);

function main_title_edit_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x_title_edit_Callback(hObject, eventdata, handles)
guidata(hObject, handles);

function x_title_edit_CreateFcn(hObject, eventdata, handles)
% This function creates the X axis title from user
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y_title_edit_Callback(hObject, eventdata, handles)
% This function creates the Y axis title from user
guidata(hObject, handles);

function y_title_edit_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

%added by Jonathan Jeffet, 2016:
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
			    set(hObject,'BackgroundColor','white');
			end
			
			
			function Residuals_y_title_edit_Callback(hObject, eventdata, handles)
			guidata(hObject, handles);
			
			% --- Executes during object creation, after setting all properties.
			function Residuals_y_title_edit_CreateFcn(hObject, eventdata, handles)
			% This function creates the Residuals' Y axis title from user
			%
			% Hint: edit controls usually have a white background on Windows.
%end of what jonathan jeffet added
			
			
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%% Change Range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function change_range_check_box_Callback(hObject, eventdata, handles)
%checkboxStatus = 0, if the box is unchecked,
%checkboxStatus = 1, if the box is checked
checkboxStatus = get(handles.change_range_check_box,'Value');
if(checkboxStatus)
    %if box is checked needs to get a range
    set(handles.x_min_edit,'Visible' , 'on');
    set(handles.x_max_edit,'Visible' , 'on');
else
    %if box is unchecked, no need for a range
    set(handles.x_min_edit,'Visible' , 'off');
    set(handles.x_max_edit,'Visible' , 'off');
end
guidata(hObject, handles)

function x_min_edit_Callback(hObject, eventdata, handles)
guidata(hObject, handles)

function x_min_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function x_max_edit_Callback(hObject, eventdata, handles)
guidata(hObject, handles)

function x_max_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function function_edit_Callback(hObject, eventdata, handles)
% Gets the name of the function to be used
function_input = get(hObject,'String');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function function_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function a0_edit_Callback(hObject, eventdata, handles)
% Gets the initial guess vector
a0 = get(hObject,'String');
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function a0_edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function help_pushbutton_Callback(hObject, eventdata, handles)
% Opens the help documentation help.doc
open('help.pdf');

%%%%%%%%%%%%%%%%%%%%%%  Saving data  %%%%%%%%%%%1%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_path_edit_Callback(hObject, eventdata, handles)
guidata(hObject, handles)

function save_path_edit_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackColor','white');
end


function save_data_pushbutton_Callback(hObject, eventdata, handles)
% saves the y_fit_plot in the required path
data = handles.Data;
x = data.x; y = data.y; dx = data.dx; dy = data.dy;

fit_data = data.Fit;
uisave({'fit_data'},'var1');
%save_path = get(handles.save_path_edit,'String'); % enter data saving path
%save(save_path,'y_fit_plot','-ascii');
guidata(hObject, handles)


% --- Executes on button press in PlotNewWindow.
function PlotNewWindow_Callback(hObject, eventdata, handles)
% hObject    handle to PlotNewWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of PlotNewWindow


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over PlotNewWindow.
function PlotNewWindow_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to PlotNewWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function PlotNewWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotNewWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function legendInput_Callback(hObject, eventdata, handles)
% hObject    handle to legendInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
legNames=get(hObject,'String');
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of legendInput as text
%        str2double(get(hObject,'String')) returns contents of legendInput as a double


% --- Executes during object creation, after setting all properties.
function legendInput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to legendInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkBoxLegend.
function checkBoxLegend_Callback(hObject, eventdata, handles)
% hObject    handle to checkBoxLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
checkboxStatusLeg = get(handles.checkBoxLegend,'Value');
if(checkboxStatusLeg)
    %if box is checked needs to get a range
    set(handles.legendInput,'Visible' , 'on');
    set(handles.pushLegend,'Visible' , 'on');
	set(handles.Legend2,'Visible' , 'on');
	set(handles.LegendText1,'Visible' , 'on');
	set(handles.LegendText2,'Visible' , 'on');
else
    %if box is unchecked, no need for a range
    set(handles.legendInput,'Visible' , 'off');
    set(handles.pushLegend,'Visible' , 'off');
	set(handles.Legend2,'Visible' , 'off');
	set(handles.LegendText1,'Visible' , 'off');
	set(handles.LegendText2,'Visible' , 'off');
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of checkBoxLegend


% --- Executes on button press in pushLegend.
function pushLegend_Callback(hObject, eventdata, handles)
% hObject    handle to pushLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%checkboxStatusLegend = get(handles.checkBoxLegend,'Value');
%if(checkboxStatusLegend)
    %if box is checked needs draw legend
    figure(1)
    legend('off')
    legend('toggle')
    [legend_h,object_h,plot_h,text_strings] =legend;

    j=1;
    legName = get(handles.legendInput,'String');
    legNames=strsplit(';', legName, 'omit');
    hVec(1) = 0;
    for i=1:length(legNames) / 2
        hVec(j) = plot_h(4 * i -1);
        hVec(j+1) = plot_h(4 * i);
        j = j +2 ;
    end
    legend(hVec , legNames);
%end
    set(legend_h, 'Interpreter', GetInterpreterString(hObject, handles, 'legendsInterpreterPopup'));


% --- Executes on button press in Legend2.
function Legend2_Callback(hObject, eventdata, handles)
% hObject    handle to Legend2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure(1)
    legend('off')
    legend('toggle')
    [legend_h,object_h,plot_h,text_strings] =legend;
    j=1;
    legName = get(handles.legendInput,'String');
    legNames=strsplit(';', legName, 'omit');
    hVec(1) = 0;
	hVec(j) = plot_h(3);
    for i=1:length(legNames)-1
        
        hVec(j+1) = plot_h(3 + i);
        j = j +1 ;
    end
    legend(hVec , legNames);
    [legend_h,object_h,plot_h,text_strings] =legend;
    set(legend_h, 'Interpreter', GetInterpreterString(hObject, handles, 'legendsInterpreterPopup'));


% --- Executes on button press in plotSavedFit.
function plotSavedFit_Callback(hObject, eventdata, handles)
% hObject    handle to plotSavedFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
{'*.mat','MATLAB Files (*.mat)';
   '*.mat','MAT-files (*.mat)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file');
 
filename_full = fullfile(pathname,filename);   
% set(handles.filename_edit,'String', filename_full);
% filename_input = get(handles.filename_edit,'String'); % gets the name of the data .dat/.txt file
data_given = load(filename_full, 'fit_data');
data_given = data_given.fit_data;
if get(handles.PlotNewWindow,'Value');
    figure()
else
    figure(1);
end
hold on
plot(data_given.x_fit_plot,data_given.y_fit_plot,'r-');        % Plots  the fitted curve
hold off


% --- Executes on button press in PasteDataButton.
function PasteDataButton_Callback(hObject, eventdata, handles)
% hObject    handle to PasteDataButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[dlg, data] = PasteDataDialog();

%if (dlg > 0) % Was the fialog successfully loaded? 
%%This statment was commented out after moving to Matlab 2015.
    uiwait(dlg);
    
    if (data.WasApproved)
        handles.Data.x = data.SelectedData(:, 1);
        handles.Data.dx = data.SelectedData(:, 2);
        handles.Data.y = data.SelectedData(:, 3);
        handles.Data.dy = data.SelectedData(:, 4);
    end
%end
1;


% --- Executes on selection change in xTitleInterpreterPopup.
function xTitleInterpreterPopup_Callback(hObject, eventdata, handles)
% hObject    handle to xTitleInterpreterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (ishandle(1))
    figure(1)
    h = get(gca, 'XLabel');
    set(h, 'Interpreter', GetInterpreterString(hObject, handles, 'xTitleInterpreterPopup'));
end


% --- Executes during object creation, after setting all properties.
function xTitleInterpreterPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xTitleInterpreterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in titleInterpreterPopup.
function titleInterpreterPopup_Callback(hObject, eventdata, handles)
% hObject    handle to titleInterpreterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (ishandle(1))
    figure(1)
    h = get(gca, 'Title');
    set(h, 'Interpreter', GetInterpreterString(hObject, handles, 'titleInterpreterPopup'));
end


% --- Executes during object creation, after setting all properties.
function titleInterpreterPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to titleInterpreterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in yTitleInterpreterPopup.
function yTitleInterpreterPopup_Callback(hObject, eventdata, handles)
% hObject    handle to yTitleInterpreterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (ishandle(1))
    figure(1)
    h = get(gca, 'YLabel');
    set(h, 'Interpreter', GetInterpreterString(hObject, handles, 'yTitleInterpreterPopup'));
end


% --- Executes during object creation, after setting all properties.
function yTitleInterpreterPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yTitleInterpreterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in legendsInterpreterPopup.
function legendsInterpreterPopup_Callback(hObject, eventdata, handles)
% hObject    handle to legendsInterpreterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (ishandle(1))
    figure(1)
    [legend_h,object_h,plot_h,text_strings] =legend;
    set(legend_h, 'Interpreter', GetInterpreterString(hObject, handles, 'legendsInterpreterPopup'));
end

% --- Executes during object creation, after setting all properties.
function legendsInterpreterPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to legendsInterpreterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% Changed by Jonathan Jeffet, 2016: adds residual y axis 
% --- Executes on selection change in ResyTitleInterperterPopup.
function ResyTitleInterperterPopup_Callback(hObject, eventdata, handles)
% hObject    handle to ResyTitleInterperterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ResyTitleInterperterPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ResyTitleInterperterPopup
if (ishandle(200))
    figure(200)
    h = get(gca, 'YLabel');
    set(h, 'Interpreter', GetInterpreterString(hObject, handles, 'ResyTitleInterperterPopup'));
end

% --- Executes during object creation, after setting all properties.
function ResyTitleInterperterPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ResyTitleInterperterPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plotChiValues.
function plotChiValues_Callback(hObject, eventdata, handles)
% hObject    handle to plotChiValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotChiValues


% --- Executes during object creation, after setting all properties.
guidata(hObject, handles)
