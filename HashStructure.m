function fHashVal = HashStructure(sStruct)

fHashVal = uint32(0);

for (nElement = 1:numel(sStruct))
   % - Get field names
   cstrFields = fieldnames(sStruct(nElement));
   
   % - Hash individual fields
   for (nField = 1:numel(cstrFields))
      % - Hash field name
      fHashVal = uint32(mod(double(fHashVal) * 10 + sum(cstrFields{nField}), intmax('uint32') - 100));
      
      % - Hash field value
      oField = sStruct(nElement).(cstrFields{nField});
      fHashVal = uint32(mod(double(fHashVal) * 10 + double(HashCellArray(CellFlatten(oField))), intmax('uint32') - 100));
   end
end

function fHashVal = HashCellArray(cArray)

% - Hash the array values
fHashVal = uint32(0);

for (nCell = 1:numel(cArray))
   % - Is this a structre?
   if (isstruct(cArray{nCell}))
      fHashVal = uint32(mod(double(fHashVal) * 10 + double(HashStructure(cArray{nCell})) * 100, intmax('uint32') - 100));
      
   elseif (isempty(cArray{nCell}))
      % - Empty array
      fHashVal = uint32(99);

   else
      % - Try to convert to number
      try
         fHashVal = uint32(mod(double(fHashVal) * 10 + sum(double(cArray{nCell})) * 100, intmax('uint32') - 100));
      catch %#ok<CTCH>
         fHashVal = uint32(99);
      end
   end
end
