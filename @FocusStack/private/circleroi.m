function [sCircle, vnIncludedPixels] = circleroi(hFigure)

% circleg - FUNCTION Define a circle on a figure
%
% Usage: [sCircle, vnIncludedPixels] = circleroi(<hFigure>)
%
% The user can define a circular ROI on a specified figure using the mouse.
%
% 'sCircle' will be a structure with fields:
%           .vfCenter   -- A vector with elements [fXCenter fYXenter], in pixels
%                          on the figure
%           .fRadius    -- A scalar containing the radius in pixels on the
%                          figure
%
% 'vnIncludedPixels' will be a vector containing linear indices into the
% currently displayed image, for pixels that fall inside the ROI.  If no image
% is displayed in the specified figure, 'vnIncludedPixels' will be empty.
%
% Uses code from roitool: http://www.mathworks.com/matlabcentral/fileexchange/9183.

% Author: Dylan Muir <dylan@ini.phys.ethz.ch>
% Created: 6th January, 2011

% - Select the specified figure
if (exist('hFigure', 'var') && ~isempty(hFigure))
   figure(hFigure);
else
   % - Otherwise select the current figure (or create one)
   figure(gcf);
end

set(gca,'units','pixels'); 

bFinishedDefinition = false;

while (~bFinishedDefinition)
   hROIHandle = get_center();
   
   % set the roisize
   disp('--- circleroi: Define the desired radius by moving the mouse, then double-click to finalise.');
   disp('        Press ''esc'' to select again.');
   
   set(gcf,'WindowButtonMotionFcn', @change_radius);
   set(gcf,'WindowButtonDownFcn', @complete_radius);
   set(gcf,'KeyPressFcn', @key_press);
   waitfor(gcf,'WindowButtonMotionFcn');
   
   if (ishandle(hROIHandle))
      userdata = get(hROIHandle, 'userdata');
   end
   drawnow;
end
   
% - Extract data 
sCircle.vfCenter = [userdata.xcenter userdata.ycenter];
sCircle.fRadius = userdata.R0;

% - Find corresponding pixels in an image, if there is one
vhImHandles = imhandles(gcf);
if (~isempty(vhImHandles))
   hImageModel = getimagemodel(vhImHandles(1));
   mbMask = poly2mask(  get(hROIHandle, 'XData'), get(hROIHandle, 'YData'), ...
                        hImageModel.getImageWidth(), hImageModel.getImageHeight());

   vnIncludedPixels = find(mbMask);
else
   vnIncludedPixels = [];
end

% - Delete line element
delete(hROIHandle);


% Sub-function - get_center
% -----------------------------------------------------------------------
function hROIHandle = get_center

   disp('--- circleroi: Select the center of the circle...');
   drawnow;
   
   [xcenter, ycenter] = ginput(1);
   x0 = xcenter; y0 = ycenter; R0 = 5;
   t = 0:pi/20:2*pi;
   xi = R0*cos(t)+xcenter;
   yi = R0*sin(t)+ycenter;
   
   hROIHandle = line(xi,yi,'LineWidth',2,'Color','red');
   userdata.t = t; userdata.xcenter = xcenter; userdata.ycenter = ycenter; userdata.R0 = R0;
   set(hROIHandle, 'userdata' ,userdata);
   
end

% Sub-function - change_radius
% -----------------------------------------------------------------------
function change_radius(src, event)  %#ok<INUSD>
   
   userdata = get(hROIHandle,'userdata');
   t = userdata.t; xcenter = userdata.xcenter; ycenter = userdata.ycenter;
   
   current_pts = get(gca,'CurrentPoint');
   current_pt = current_pts(1,1:2);
   R0 = sqrt((current_pt(1)-xcenter)^2+(current_pt(2)-ycenter)^2);
   xi = R0*cos(t)+ xcenter;
   yi = R0*sin(t)+ ycenter;
   userdata.R0 = R0;
   set(hROIHandle,'Xdata',xi,'Ydata',yi);
   set(hROIHandle,'userdata',userdata);
   drawnow;
   mouseclick = get(gcf,'SelectionType');

end

% Sub-function - complete_radius
% -----------------------------------------------------------------------
function complete_radius(src, event) %#ok<INUSD>
   mouseclick = get(gcf,'SelectionType');
   if strcmp(mouseclick,'open')
      bFinishedDefinition = true;
      set(gcf,'WindowButtonMotionFcn','');
      set(gcf,'WindowButtonDownFcn','');
   end
end

% Sub-fucntion - key_press
% -----------------------------------------------------------------------
function key_press(src, event) %#ok<INUSL>
   
   switch event.Key,  %process shortcut keys
      case 'escape'
         % - Cancel radius selection
         set(gcf,'WindowButtonMotionFcn', []);
         
         % - Delete ROI marker
         delete(hROIHandle);
         
      otherwise
         % - Do nothing
   end
end

end
   
% --- END of circleg.m ---
