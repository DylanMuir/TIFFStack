function [sRegions] = ROIs2Regions(cvsROIs, vnImageSize)

% - Build a structure
sRegions.Connectivity = 8;
sRegions.ImageSize = vnImageSize;
sRegions.NumObjects = numel(cvsROIs);
sRegions.PixelIdxList = {};

for (nROIIndex = numel(cvsROIs):-1:1)
   sThisROI = cvsROIs{nROIIndex};
   
   switch (lower(sThisROI.strType))
      case 'rectangle'
         if (isfield(sThisROI, 'strSubtype') && isequal(lower(sThisROI.strSubtype), 'shape'))
            % - Skip this one
            break;
            
         else
            % - Make a rectangular mask
            mbThisMask = false(vnImageSize);
            mbThisMask(sThisROI.vnRectBounds(2):sThisROI.vnRectBounds(4), sThisROI.vnRectBounds(1):sThisROI.vnRectBounds(3)) = true;
            Regions.PixelIdxList{nROIIndex} = find(mbThisMask);
         end
         
      case 'oval'
         % - Draw an oval inside the bounding box
         mbThisMask = ellipse2mask('bounds', vnImageSize, sThisROI.vnRectBounds);
         sRegions.PixelIdxList{nROIIndex} = find(mbThisMask);
         
      otherwise
         warning( 'ROIs2Regions:unsupported', ...
                  '--- ROIs2Regions: Warning: Unsupported ROI type.');
   end
end

