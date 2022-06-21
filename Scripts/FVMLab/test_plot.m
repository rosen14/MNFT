## Useful since Octave 4.0

close all
clear h

graphics_toolkit qt

h.ax = axes ("position", [0.05 0.42 0.5 0.5]);

function update_plot (obj, init = false)

  ## gcbo holds the handle of the control
  h = guidata (obj);

  a = get (h.slider1, "value");
  x_range = [-2:0.2:2];
  y_range = [-2:0.2:2];
  [x,y] = meshgrid(x_range, y_range);
  frames = 500;
  z = a*(x) .* y;

  if (init)
    h.plot = surf(x, y, z);
  else
    set (h.plot, 'ZData', z);
    set (h.plot, 'CData', z);
  endif
  guidata (obj, h);

endfunction

h.slider1 = uicontrol ("style", "slider",
                      "units", "normalized",
                      "string", "slider",
                      "callback", @update_plot,
                      "value", 0.4,
                      "position", [0.05 0.25 0.35 0.06]);



set (gcf, "color", get(0, "defaultuicontrolbackgroundcolor"))
guidata (gcf, h)
update_plot (gcf, true);
