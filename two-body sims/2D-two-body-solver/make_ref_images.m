function make_ref_images(csv_dir, csv_filename, tag, x0, tf, params)
    images_dir = [csv_dir '_IMAGES'];
    plottable_dir = [csv_dir '_' tag];

    % Get OSI
    A = readmatrix([csv_dir,'/',csv_filename]);
    A = sortrows(A);
    path_struct.mu = A(:,1)';
    path_struct.p1 = A(:,2)';
    path_struct.p2 = A(:,3)';
    path_struct.p3 = A(:,4)'; 
    path_struct.p4 = A(:,5)';
    path_struct.p5 = A(:,6)';
    path_struct.p6 = A(:,7)';
    path_struct.num_conj_pt = A(:,8)';

    % Find OSI
    num_conj_diff = path_struct.num_conj_pt(2:end) - path_struct.num_conj_pt(1:(end-1));
    osi_body = find(num_conj_diff ~= 0) + 1;
    path_struct.osi = [1 osi_body length(path_struct.mu)];

    header = {};

    % Get the ivp solution for each OSI pair (denoting a segment of the total path)
    for i = 1:length(path_struct.osi)-1
        if path_struct.mu(path_struct.osi(i)) < 0.001
            % On the first iteration, show the image for mu=0, not the 
            % halfway point, drawn a quarter of the way through the segment
            ivp_idx = path_struct.osi(i);
            quarter_seg = floor((path_struct.osi(i+1) + path_struct.osi(i)) / 4);
            if quarter_seg == 0
                quarter_seg= 1;
            end
            render_coords = [path_struct.p1(quarter_seg) path_struct.p2(quarter_seg) ...
                path_struct.p3(quarter_seg) path_struct.p4(quarter_seg) ...
                path_struct.p5(quarter_seg) path_struct.p6(quarter_seg)];
        else
            if i ~= 1 && i ~= length(path_struct.osi)-1 && ...
                    i ~= ceil((length(path_struct.osi)-1)/2)
                continue
            end
            ivp_idx = floor((path_struct.osi(i+1) + path_struct.osi(i)) / 2);
            render_coords = [path_struct.p1(ivp_idx) path_struct.p2(ivp_idx) ...
                path_struct.p3(ivp_idx) path_struct.p4(ivp_idx) ...
                path_struct.p5(ivp_idx) path_struct.p6(ivp_idx)];
        end
        p0 = [path_struct.p1(ivp_idx) path_struct.p2(ivp_idx) ...
            path_struct.p3(ivp_idx) path_struct.p4(ivp_idx) ...
            path_struct.p5(ivp_idx) path_struct.p6(ivp_idx)];
        ivp_solution = solve_IVP(x0,p0,tf,params,path_struct.mu(ivp_idx),1);

        % Plot that solution
%        fprintf('plotting p0=%s for mu=%f...\n', ...
%            sprintf('%.2f ',ivp_solution.p(1,:)),path_struct.mu(ivp_idx));

        figure('visible', 'off');
        plot(ivp_solution.x(:,1),ivp_solution.x(:,2),'b-','LineWidth',5); % Solution curve 1
        hold on
        plot(ivp_solution.x(:,4),ivp_solution.x(:,5),'c-','LineWidth',5); % Solution curve 2
        plot(ivp_solution.x(end,1),ivp_solution.x(end,2),'ro'); % Endpoint 1
        plot(ivp_solution.x(end,4),ivp_solution.x(end,5),'ro'); % Endpoint 2
        plot(ivp_solution.x(1,1),ivp_solution.x(1,2),'go'); % Start point 1
        plot(ivp_solution.x(1,4),ivp_solution.x(1,5),'go'); % Start point 2
        hold off
        ax = gca;
        set(ax, 'Color', 'none');
        set(ax,'Visible','off');
        daspect([1 1 1]);

        % Save solution as a png
        image_filename = [csv_filename(1:(end-4)), ...
            sprintf('(seg=%.2f).png',path_struct.mu(ivp_idx))];
        %print([images_directory,'/',image_filename],'-dpng');
        
        cdir= cd('altmany-export_fig-b1a7288');
        export_fig(['../',images_dir,'/',image_filename], '-transparent');
        cd(cdir);

        plot_function(ivp_solution, ivp_solution.x(end,:));

        % Generate header text
        image_ref = strcat(images_dir,'/',image_filename, ...
            sprintf(',%f',render_coords));
        header(end+1) = {image_ref};
    end
    
    % Rewrite body text with header back to file
    body = fileread([csv_dir,'/',csv_filename]);
%    fprintf('writing to %s...\n', [plottable_dir,'/','header_',csv_filename]);
    fd = fopen([plottable_dir,'/','header_',csv_filename],'w');
    fprintf(fd,'%d\n',length(header));
    fprintf(fd,'%s\n',header{:});
    fprintf(fd,'%s\n',body);
    fclose(fd);
%    fprintf('done!\n');

end