function [fDCOffset, fPhase, fFrequency, fAmplitude] = fit_sin(vfAngles, vfResponses, bDisplay)

   if (~exist('bDisplay', 'var'))
      bDisplay = false;
   end

   if (~any(isnan(vfResponses)))
      % - Initial params
      fDCOffset = mean(vfResponses);
      fFrequency = 4;
      fPhase = pi;
      fAmplitude = max(vfResponses)-min(vfResponses);
      
      x0 = [fDCOffset fPhase fFrequency fAmplitude];
      if (bDisplay)
         vfVisAngles = linspace(min(vfAngles), max(vfAngles), 200);
      end
      
      vfInterpAngles = linspace(min(vfAngles), max(vfAngles), 32);
      vfInterpResponses = interp1(vfAngles, vfResponses, vfInterpAngles);
      
      x = fminunc(@sin_cost,x0);
   
   else
      x = nan;
   end
   
   % - Return zero vectors for nan results
   if (any(isnan(x)))
      x = zeros(4, 1);
   end
   
   % - Extract parameters
   fDCOffset = x(1);
   fPhase = x(2);
   fFrequency = x(3);
   fAmplitude = x(4);
   

   function [fSSErr] = sin_cost(x)
      % - Extract parameters
      fDCOffset = x(1);
      fPhase = x(2);
      fFrequency = x(3);
      fAmplitude = x(4);
      
      % - Compute sin curve for supplied parameters
      vfFun = fAmplitude * sin(vfInterpAngles .* fFrequency + fPhase) + fDCOffset;
      
      if (bDisplay)
         vfVisFun = fAmplitude * sin(vfVisAngles .* fFrequency + fPhase) + fDCOffset;
         
         % - Visualise fitting process
         figure(101);
         clf;
         plot(vfInterpAngles, vfInterpResponses, 'k-');
         hold on;
         plot(vfVisAngles, vfVisFun, 'b--');
         drawnow;
      end
      
      % - Compute sum-squared error
      fSSErr = sum((vfInterpResponses - vfFun).^2);
   end

end
