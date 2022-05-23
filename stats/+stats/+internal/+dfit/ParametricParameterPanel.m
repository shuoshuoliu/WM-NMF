 classdef ParametricParameterPanel < stats.internal.dfit.ParameterPanel
    % ParametricParmeterPanel 
    
    % Copyright 2019-2020 The MathWorks, Inc

    properties (Access = private)
       NumParams;
       FixedValueEditFields;
       UseEstimated;
       FitObj;
   end
   
   methods
       
       function this = ParametricParameterPanel(fo, panelGrid, gridRow, gridColumn, ft)          
           this.FitObj = fo;
           this.FitType = ft;
           this.DisplayName = ft.getDisplayName();
           
           this.Panel = uipanel(panelGrid, 'Title', this.DisplayName, 'Visible', 'off');
           this.Panel.Layout.Row = gridRow;
           this.Panel.Layout.Column = gridColumn;
                  
           pnames = ft.getParameterNames();
           pdesc = ft.getParameterDescriptions();
           preq = ft.getParameterRequired();
                
           g = uigridlayout(this.Panel, [length(pnames) + 1, 2]);
           g.ColumnWidth = {'fit', '1x'};
           rowCell = cell(1, length(pnames) + 1);
           for i = 1:(length(pnames) + 1)
               rowCell{1, i} = 'fit';
           end
           g.RowHeight = rowCell;
           
           l = uilabel(g, 'Text', [getString(message('stats:dfittool:label_distributionParameters')) ' ']);
           l.Layout.Row = 1;
           l.Layout.Column = [1 2];
           this.NumParams = length(pnames);
               
           this.FixedValueEditFields = cell(1, this.NumParams);
           this.UseEstimated = zeros(1, this.NumParams);
           for i = 1:this.NumParams
               if preq(i)
                   l = uilabel(g, 'Text', ['  ' pnames{i}, ' (' pdesc{i} ') ' getString(message('stats:dfittool:label_requiresFixed'))]);
                   l.Layout.Row = i + 1;
                   l.Layout.Column = 1;
                   this.FixedValueEditFields{i} = uieditfield(g, 'ValueChangingFcn', @(~, ~)this.cbEFChanging());
                   this.FixedValueEditFields{i}.Layout.Row = i+1;
                   this.FixedValueEditFields{i}.Layout.Column = 2;
                   this.UseEstimated(i) = false;
               else
                   this.UseEstimated(i) = true;
                   l = uilabel(g, 'Text', ['  ' pnames{i} ' (' pdesc{i} ')']);
                   l.Layout.Row = i + 1;
                   l.Layout.Column = [1 2];
               end
           end
       end
       
       function tf = isBadValue(this)
           err = '';
           pnames = this.FitType.getParameterNames();
           preq = this.FitType.getParameterRequired();
           for i = 1:this.NumParams
               if preq(i)
                   value = this.FixedValueEditFields{i}.Value;
                   if isnan(str2double(value)) % error
                       err = [err getString(message('stats:dfittool:error_badFixedValue', pnames{i}))];  %#ok<AGROW>
                   end
               end
           end
           if isempty(err)
               tf = false;
           else
               tf = true;
               stats.internal.dfit.errordlg(err, getString(message('stats:dfittool:title_invalidValue')));
           end
       end
       
	   function updateFixedValues(this, newUDDObj)
           for i = 1:length(this.FixedValueEditFields)
               if isa(this.FixedValueEditFields{i}, 'matlab.ui.control.EditField')
                    this.FixedValueEditFields{i}.Value =  newUDDObj.pfixedtext{i};  
               end
           end
       end
       
       function fitArgs = getFitArgs(pp, fitType, fitEditor)
           fitArgs = {fitType.getModelType, getFit(fitEditor), getFitName(fitEditor), fitType.getRealName(), getDataName(fitEditor), fitEditor, fitEditor.getExclusionRule(), pp.UseEstimated, pp.getFixedVals()}; 
       end
              
   end
   
   methods(Access = private)
       function cbEFChanging(this)
            this.FitObj.clearResults;
            this.FitObj.setFitChanged();
       end

       function fixedVals = getFixedVals(this)
           fixedVals = cell(1, this.NumParams);
           for i = 1:this.NumParams
               if isa(this.FixedValueEditFields{i}, 'matlab.ui.control.EditField')
                   fixedVals{i} = this.FixedValueEditFields{i}.Value;
               else
                   fixedVals{i} = '';
               end
           end
       end
   end
end
