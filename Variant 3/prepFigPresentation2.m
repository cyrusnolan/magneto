function prepFigPresentation2(fignum)
        %
        % prepares a figure for presentations
        %
        % Fontsize: 14
        % Fontweight: bold
        % LineWidth: 1.5
        
            figure(fignum);
            fig_children=get(fignum,'children'); %find all sub-plots
            
            for i=1:length(fig_children)
            
                if class(fig_children(i)) == "matlab.graphics.illustration.Legend"
                    set(fig_children(i),'FontSize',16,'Interpreter','latex');
                end

                set(fig_children(i),'FontSize',16);
                % set(fig_children(i),'FontWeight','bold');
                
                fig_children_children=get(fig_children(i),'Children');
                set(fig_children_children,'LineWidth',1.5);
            end
        end