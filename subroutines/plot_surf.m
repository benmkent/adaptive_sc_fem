function plot_surf(ax,u, fem, params)

switch params.grid
    case 'q1'
        [u_reshape, x,y] = vec2xy(u, fem);
        cla(ax);
        switch params.plottype
            case 'surf'
                surf(ax,x, y, (u_reshape));
                view(ax,params.plotview);
            case 'contourf'
                contourf(ax,x, y, u_reshape,20);
        end
        colorbar(ax);
    case 'p1'
        cla(ax);
        axes(ax);
        switch  params.plottype
            case 'surf'
                trisurf(fem.T, fem.xy(:,1), fem.xy(:,2), abs(u));
                view(ax,params.plotview);
            case 'contourf'
                [u_reshape, x,y] = vec2xy(u, fem);
                contourf(ax,x, y, (u_reshape));
        end
        colorbar(ax);
end

end