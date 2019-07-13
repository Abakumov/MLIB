function xx=x2grid(x, Gox, Gdx, Gnx)
xx = round((x-Gox)/Gdx)+1;
xx = max(1, xx);
xx = min(xx, Gnx);
