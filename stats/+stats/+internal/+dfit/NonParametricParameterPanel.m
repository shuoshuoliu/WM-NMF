classdef NonParametricParameterPanel < stats.internal.dfit.ParameterPanel
    % NonParametricParmeterPanel 
    
    % Copyright 2019-2020 The MathWorks, Inc
    
   properties (Access = private)
        KernelDropdown
        SpecifyBandWidthEditfield
        SpecifyDomainLowBoundEditfield
        SpecifyDomainHighBoundEditfield
        BandwidthAutoRadioButton
        BandwidthSpecifyRadioButton
        DomainUnboundedRadioButton
        DomainPositiveRadioButton
        DomainSpecifyRadioButton
        FitObj
   end
   
   properties (Constant)
        BANDWIDTH_AUTO = 0;
        BANDWIDTH_SPECIFY = 1;

        DOMAIN_UNBOUNDED = 0;
        DOMAIN_POSITIVE = 1;
        DOMAIN_SPECIFY = 2;
   end
                   
   methods
       
       function this = NonParametricParameterPanel(fo, panelGrid, gridRow, gridColumn, ft)            
            this.FitObj = fo;
            this.FitType = ft;
            this.DisplayName = ft.getDisplayName();
            
            this.Panel = uipanel(panelGrid, 'Title', this.DisplayName, 'Visible', 'off');
            this.Panel.Layout.Row = gridRow;
            this.Panel.Layout.Column = gridColumn;
            mainNPPGrid = uigridlayout(this.Panel, [6, 3]);
            mainNPPGrid.RowHeight = {'fit', 18, 18, 18, 18, 18};
            mainNPPGrid.ColumnWidth = {'fit', 24, '1x'};
            
            l = uilabel(mainNPPGrid, 'Text', getString(message('stats:dfittool:label_kernel')));
            l.Layout.Row = 1;
            l.Layout.Column = 1;
            
            l = uilabel(mainNPPGrid, 'Text', getString(message('stats:dfittool:label_bandwidth')));
            l.Layout.Row = 2;
            l.Layout.Column = 1;
            
            l = uilabel(mainNPPGrid, 'Text', getString(message('stats:dfittool:label_domain')));
            l.Layout.Row = 4;
            l.Layout.Column = 1;              
        
            this.KernelDropdown = uidropdown(mainNPPGrid, 'Items',...
                {getString(message('stats:dfittool:kernelType_normal')), ...
                getString(message('stats:dfittool:kernelType_box')), ...
                getString(message('stats:dfittool:kernelType_triangle')), ...
                getString(message('stats:dfittool:kernelType_epanechnikov'))}, ...
                'ValueChangedFcn', @(~, ~)this.cbButtonOrDropdownChanged());
            this.KernelDropdown.Layout.Row = 1;
            this.KernelDropdown.Layout.Column = [2, 3];
            
            l = uilabel(mainNPPGrid, 'Text', getString(message('stats:dfittool:radiobutton_auto')));
            l.Layout.Row = 2;
            l.Layout.Column = 3;
            
            specifyBandWidthGrid = uigridlayout(mainNPPGrid, [1, 2]);
            specifyBandWidthGrid.Padding = [0 0 0 0];
            specifyBandWidthGrid.RowHeight = {18};
            specifyBandWidthGrid.ColumnWidth = {'fit', '1x'};
            specifyBandWidthGrid.Layout.Row = 3;
            specifyBandWidthGrid.Layout.Column = 3;
            
            l = uilabel(specifyBandWidthGrid, 'Text', getString(message('stats:dfittool:radiobutton_specify')));
            l.Layout.Row = 1;
            l.Layout.Column = 1;
            
            this.SpecifyBandWidthEditfield = uieditfield(specifyBandWidthGrid, 'ValueChangingFcn', @(~, ~)this.cbBandWidthEditfieldChanging());
            this.SpecifyBandWidthEditfield.Layout.Row = 1;
            this.SpecifyBandWidthEditfield.Layout.Column = 2;
            
            l = uilabel(mainNPPGrid, 'Text', getString(message('stats:dfittool:radiobutton_unbounded')));
            
            l.Layout.Row = 4;
            l.Layout.Column = 3;
            
            l = uilabel(mainNPPGrid, 'Text', getString(message('stats:dfittool:radiobutton_positive')));
            l.Layout.Row = 5;
            l.Layout.Column = 3;
            
            specifyDomainGrid = uigridlayout(mainNPPGrid, [1, 4]);
            specifyDomainGrid.Padding = [0 0 0 0];
            specifyDomainGrid.RowHeight = {18};
            specifyDomainGrid.ColumnWidth = {'fit', '1x', 'fit', '1x'};
            specifyDomainGrid.Layout.Row = 6;
            specifyDomainGrid.Layout.Column = 3;
            
            l = uilabel(specifyDomainGrid, 'Text', getString(message('stats:dfittool:radiobutton_specify')));
            l.Layout.Row = 1;
            l.Layout.Column = 1;
            
            this.SpecifyDomainLowBoundEditfield = uieditfield(specifyDomainGrid, 'ValueChangingFcn', @(~, ~)this.cbDomainEditfieldChanging());
            this.SpecifyDomainLowBoundEditfield.Layout.Row = 1;
            this.SpecifyDomainLowBoundEditfield.Layout.Column = 2;
            
            l = uilabel(specifyDomainGrid, 'Text', [' ' getString(message('stats:dfittool:label_to')) ' ']);
            l.Layout.Row = 1;
            l.Layout.Column = 3;
            
            this.SpecifyDomainHighBoundEditfield = uieditfield(specifyDomainGrid, 'ValueChangingFcn', @(~, ~)this.cbDomainEditfieldChanging());
            this.SpecifyDomainHighBoundEditfield.Layout.Row = 1;
            this.SpecifyDomainHighBoundEditfield.Layout.Column = 4;
            
            bg1 = uibuttongroup(mainNPPGrid, 'BorderType', 'none', 'SelectionChangedFcn', @(~, ~)this.cbButtonOrDropdownChanged());
            bg1.Layout.Row = [2 3];
            bg1.Layout.Column = 2;
            
            bg2 = uibuttongroup(mainNPPGrid, 'BorderType', 'none', 'SelectionChangedFcn', @(~, ~)this.cbButtonOrDropdownChanged());
            bg2.Layout.Row = [4 6];
            bg2.Layout.Column = 2;
            
            this.BandwidthAutoRadioButton = uiradiobutton(bg1, 'Text', '', 'Position', [10 28 15 20]);
            this.BandwidthSpecifyRadioButton = uiradiobutton(bg1, 'Text', '', 'Position', [10 0 15 20]);
            this.DomainUnboundedRadioButton = uiradiobutton(bg2, 'Text', '', 'Position', [10 56 15 20]);
            this.DomainPositiveRadioButton = uiradiobutton(bg2, 'Text', '', 'Position', [10 28 15 20]);
            this.DomainSpecifyRadioButton = uiradiobutton(bg2, 'Text', '', 'Position', [10 0 15 20]);
       end
       
       function tf = isBadValue(this)
           err = '';

           if (getBandwidthType(this) == this.BANDWIDTH_SPECIFY)
               value = this.SpecifyBandWidthEditfield.Value;
               if isnan(str2double(value)) % error
                   err = [err getString(message('stats:dfittool:error_badBandWidth'))];  
               end
           end
           if (getDomainType(this) == this.DOMAIN_SPECIFY)
               validLowBounds = true;
               validHighBounds = true;
               lowBounds = -Inf;
               highBounds = Inf;
               
               value = this.SpecifyDomainLowBoundEditfield.Value;
               if isnan(str2double(value)) % error
                   validLowBounds = false;
                   err = [err getString(message('stats:dfittool:error_badLowBound'))]; 
               else
                   lowBounds = str2double(value);
               end
               value = this.SpecifyDomainHighBoundEditfield.Value;
               if isnan(str2double(value)) % error
                   validHighBounds = false;
                   err = [err getString(message('stats:dfittool:error_badHighBound'))]; 
               else
                   highBounds = str2double(value);
               end
               
               if (validLowBounds && validHighBounds && lowBounds >= highBounds)
                   err = [err getString(message('stats:dfittool:error_lowBndsMustBeLessThanHighBnds'))];
               end
           end

           if isempty(err)
               tf = false;
           else
               tf = true;
               stats.internal.dfit.errordlg(err, getString(message('stats:dfittool:title_invalidValue')));
           end
       end
       
       function fitArgs = getFitArgs(pp, fitType, fitEditor)
            kernelType = iGetKernelTypeFromKernelTypeDisplayName(pp.KernelDropdown.Value);
            fitArgs = {fitType.getModelType, getFit(fitEditor), getFitName(fitEditor), ...
                kernelType, pp.getBandwidthType(), ...
                pp.SpecifyBandWidthEditfield.Value, getDataName(fitEditor), fitEditor, pp.getDomainType() ...
                pp.getDomainBounds(), fitEditor.getExclusionRule()}; 
       end
       
       function updateFixedValues(this, newUddObject)
           this.KernelDropdown.Value = iGetKernelTypeDisplayNameFromKernelType(newUddObject.kernel);
           if (newUddObject.bandwidthradio == 0)
                this.BandwidthAutoRadioButton.Value = true;
           else
                this.BandwidthSpecifyRadioButton.Value = true;
           end
           this.SpecifyBandWidthEditfield.Value = newUddObject.bandwidthtext;
           if (newUddObject.supportradio == 0)
               this.DomainUnboundedRadioButton.Value = true;
           elseif (newUddObject.supportradio == 1)
               this.DomainPositiveRadioButton.Value = true;
           else
               this.DomainSpecifyRadioButton.Value = true;
           end
           this.SpecifyDomainLowBoundEditfield.Value = newUddObject.supportlower;
           this.SpecifyDomainHighBoundEditfield.Value = newUddObject.supportupper;
       end
   end
   methods(Access = private)

       function bandwithType = getBandwidthType(this)
           if (this.BandwidthAutoRadioButton.Value == 1)
               bandwithType = this.BANDWIDTH_AUTO;
           else
               bandwithType = this.BANDWIDTH_SPECIFY;
           end
       end
              
       function domainType = getDomainType(this)
            if (this.DomainUnboundedRadioButton.Value == 1)
                domainType = this.DOMAIN_UNBOUNDED;
            elseif (this.DomainPositiveRadioButton.Value == 1)
                domainType =  this.DOMAIN_POSITIVE;
            else
                domainType = this.DOMAIN_SPECIFY;
            end
       end
             
       function domainBounds = getDomainBounds(this)
           domainBounds{1} = this.SpecifyDomainLowBoundEditfield.Value;
           domainBounds{2} = this.SpecifyDomainHighBoundEditfield.Value;
       end
       
       function cbDomainEditfieldChanging(this)
            this.FitObj.clearResults;
            this.FitObj.setFitChanged();
            this.DomainSpecifyRadioButton.Value = 1;
       end
       
       function cbBandWidthEditfieldChanging(this)
            this.FitObj.clearResults;
            this.FitObj.setFitChanged();
            this.BandwidthSpecifyRadioButton.Value = 1;
       end
       
       function cbButtonOrDropdownChanged(this)
            this.FitObj.clearResults;
            this.FitObj.setFitChanged();
       end
   end
end

% helper functions
function kernelType = iGetKernelTypeFromKernelTypeDisplayName(kernalTypeDisplayName)
    switch kernalTypeDisplayName
        case getString(message('stats:dfittool:kernelType_normal'))
            kernelType = 'normal';
        case getString(message('stats:dfittool:kernelType_box'))
            kernelType = 'box';
        case getString(message('stats:dfittool:kernelType_triangle'))
            kernelType = 'triangle';
        case getString(message('stats:dfittool:kernelType_epanechnikov'))
            kernelType = 'epanechnikov';
        otherwise
            kernelType = 'normal';
    end
end

function kernelTypeDisplayName = iGetKernelTypeDisplayNameFromKernelType(kernalType)
    switch kernalType
        case 'normal'
            kernelTypeDisplayName =  getString(message('stats:dfittool:kernelType_normal'));
        case 'box'
            kernelTypeDisplayName =  getString(message('stats:dfittool:kernelType_box'));
        case 'triangle'
            kernelTypeDisplayName =  getString(message('stats:dfittool:kernelType_triangle'));
        case 'epanechnikov'
            kernelTypeDisplayName =  getString(message('stats:dfittool:kernelType_epanechnikov'));
        otherwise
            kernelTypeDisplayName = getString(message('stats:dfittool:kernelType_normal'));
    end
end