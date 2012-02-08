function [mbMask] = ellipse2mask(strType, vnImageSize, varargin)

% ellipse2mask - FUNCTION Make an ellipse mask
%
% Usage: [mbMask] = ellipse2mask('bounds',    vnImageSize, [fX1 fY1 fX2 fY2])
%                   ellipse2mask('majoraxis', vnImageSize, [fX1 fY1 fX2 fY2], fAspectRatio)
%                   ellipse2mask('center',    vnImageSize, [fCx fCy], [fMajorAxisLength fMinorAxisLength], fThetaRad)
%                   ellipse2mask('axes',      vnImageSize, Major:[fX1 fY1 fX2 fY2], Minor:[fX1 fY1 fX2 fY2])
%
% 'vnImageSize' is the size of the desired output mask, with the format [nRows
% nColumns].
%
% 'bounds': Define an ellipse within a rectangle aligned with the X and Y axis.
% The bounds of the rectangle are defined with the vector [fX1 fY1 fX2 fY2].
%
% 'majoraxis': Define an ellipse to fit within a defined major axis, along the
% line defined by [fX1 fY1 fX2 fY2], with a specified major:minor axis length
% aspect ratio 'fAspectRatio'.
%
% 'center': Define an ellipse by providing a center [fCx fCy], a major and minor
% axis length [fMajorAxisLength fMinorAxisLength] and a rotation of the major
% axis from the X axis, in radians.
%
% 'axes': Define an ellipse by providing the major and minor axis lines,
% Major:[fX1 fY1 fX2 fY2], Minor:[fX1 fY1 fX2 fY2].  WARNING: ellipse2mask will
% not fix your incorrectly defined axes!  The major and minor axes should
% intersect at their middles, and be orthogonal.

% Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 10th August, 2011

% -- Check arguments

if (nargin == 0)
   disp('*** ellipse2mask: Incorrect usage.');
   help ellipse2mask;
   return;
end

% - Is the image size a 2-element vector?
if (numel(vnImageSize) ~= 2)
   error('ellipse2mask:arguments', ...
      '''vnImageSize'' must have two elements.');
end

% - Was the function call reasonable?
bIncorrectCall = false;
switch lower(strType)
   case 'bounds'
      if (nargin < 3)
         bIncorrectCall = true;
      end
      
      % - Extract arguments and check
      vfRectBounds = varargin{1};
      
      if (numel(vfRectBounds) ~= 4)
         bIncorrectCall = true;
      end
            
   case 'majoraxis'
      if (nargin < 4)
         bIncorrectCall = true;
      end
      
      % - Extract arguments and check
      vfMajorAxis = varargin{1};
      fAspectRatio = varargin{2};

      if ((numel(vfMajorAxis) ~= 4) || ~isscalar(fAspectRatio))
         bIncorrectCall = true;
      end
      
   case 'center'
      if (nargin < 5)
         bIncorrectCall = true;
      end
      
      % - Extract arguments and check
      vfCenter = varargin{1};
      vfAxes = varargin{2};
      fTheta = varargin{3};
      
      if ((numel(vfCenter) ~= 2) || (numel(vfAxes) ~= 2) || ~isscalar(fTheta))
         bIncorrectCall = true;
      end
  
   
   case 'axes'
      if (nargin < 4)
         bIncorrectCall = true;
      end
      
      % - Extract arguments and check
      vfAAxisLine = varargin{1};
      vfBAxisLine = varargin{2};

      if ((numel(vfAAxisLine) ~= 4) || (numel(vfBAxisLine) ~= 4))
         bIncorrectCall = true;
      end
   
   otherwise
      error('ellipse2mask:unknown_type', ...
            '*** ellipse2mask: Unknown ellipse definition type.');
end

% - No, so generate an error
if (bIncorrectCall)
   error('ellipse2mask:arguments', ...
         '*** ellipse2mask: Incorrect usage.');
end


% -- Convert parameters to center plus semi-axes vectors and lengths
% vfCenter: [fCX fCY]
% vfAAxis: Unit vector along major axis
% cfBAxis: Unit vector along minor axis
% fA: Semi-major axis length
% fB: Semi-minor axis length

switch lower(strType)
   case 'bounds'
      % - Make an ellipse that fits inside an axis-aligned rectangle
      vfCenter = [vfRectBounds(1) + vfRectBounds(3) vfRectBounds(2) + vfRectBounds(4)] / 2;
      vfAAxis = [1 0];
      % - Make fA and fB one unit bigger, to fit nicely in bounds
      fA = abs(vfRectBounds(3) - vfCenter(1)) + 1;
      fB = abs(vfRectBounds(4) - vfCenter(2)) + 1;
      
   case 'majoraxis'
      % - Make an ellipse that fits along a defined major axis, with a specified
      %     major:minor axis ratio
      vfCenter = [vfMajorAxis(1) + vfMajorAxis(3) vfMajorAxis(2) + vfMajorAxis(4)] / 2;
      vfAAxis =  vfMajorAxis(3:4) - vfCenter;
      fA = sqrt(sum(vfAAxis.^2));
      vfAAxis = vfAAxis ./ fA;
      fB = fA ./ fAspectRatio;
      
   case 'center'
      % - Make an ellipse defined by a center point, two axis lengths (major and
      %     minor) and a rotation
      fA = vfAxes(1)/2;
      fB = vfAxes(2)/2;
      vfAAxis = [cos(fTheta) sin(fTheta)];
      
   case 'axes'
      % - Make an ellipse defined by the major and minor axis lines.
      vfCenter = [vfAAxisLine(1) + vfAAxisLine(3) vfAAxisLine(2) + vfAAxisLine(4)] / 2;
      vfCenterB = [vfBAxisLine(1) + vfBAxisLine(3) vfBAxisLine(2) + vfBAxisLine(4)] / 2;
      vfAAxis = vfAAxisLine(3:4) - vfAAxisLine(1:2);
      vfBAxis = vfBAxisLine(3:4) - vfBAxisLine(1:2);
      fA = sqrt(sum(vfAAxis.^2))/2;
      fB = sqrt(sum(vfBAxis.^2))/2;
      vfAAxis = vfAAxis ./ fA;
      vfBAxis = vfBAxis ./ fB;
      
      % - Check to see if the axis definitions match
      if (sum(vfBAxis - [-vfAAxis(2) vfAAxis(1)]) > 2*eps)
         warning( 'ellipse2mask:axismismatch', ...
                  '--- ellipse2mask: Warning: The major and minor axes are not orthogonal.  The output is undetermined.');
      end
      
      if (sum(vfCenter - vfCenterB) > 2*eps)
         warning( 'ellipse2mask:axismismatch', ...
                  '--- ellipse2mask: Warning: The major and minor axes do not intersect at their middles.  The output is undetermined.');
      end
end


% -- Find points within the ellipse

% - Make a mesh of points
[mnY, mnX] = meshgrid(1:vnImageSize(1), 1:vnImageSize(2));

% - Transform points
mnXT = (mnX - vfCenter(1)).*vfAAxis(1) + (mnY - vfCenter(2)).*vfAAxis(2);
mnYT = -(mnX - vfCenter(1)).*vfAAxis(2) + (mnY - vfCenter(2)).*vfAAxis(1);

% - Find points inside the ellipse
mbMask = (mnXT.^2 / fA^2 + mnYT.^2 / fB^2) < 1;

% --- END of ellipse2mask.m ---
