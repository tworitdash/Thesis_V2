R_axis = [R(1) R];

for i = 1:n
    L_axis(i) = sum(L(1:i));
end
L_axis = [0 L_axis];
figure;
plot(L_axis, R_axis, 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840] );
hold on;
plot(L_axis, -R_axis, 'LineWidth', 2, 'Color', [0.6350, 0.0780, 0.1840] );
title(['Optimum Horn, Aprture Radius = ', num2str(R(end)./lamb), ' \lambda ', 'Horn Length = ', num2str(Len./lamb), ' \lambda'], 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Horn Length cut [m]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Horn Radius cut [m]', 'FontSize', 12, 'FontWeight', 'bold');
grid on;