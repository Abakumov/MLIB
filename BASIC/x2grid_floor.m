function xx=x2grid_floor(x, Gox, Gdx, Gnx)

xx = floor((x-Gox)/Gdx)+1;
xx = max(1, xx);
xx = min(xx, Gnx);