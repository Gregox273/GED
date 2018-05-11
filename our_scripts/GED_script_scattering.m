% pick material index
material_index = 5;

% figure setup
clf
figure(1)

% get energy and mu data and plot
E = material.mev;
mu = material.coeffs(:,material_index);
plot(E, mu)

% now for log scales
figure(2)
plot(log(E),log(mu))
hold on

% get estimate of y intercept for photoelectric scattering
log_K_photo = mean(log(mu(1:10)) + 3 * log(E(1:10)))

x = log(material.mev);
y = -3*log(material.mev)+log_K_photo;
plot(x,y)
scatter(x(10), y(10))

% get estimate of y intercept for compton scattering
log_K_compton = mean(log(mu(50:200)) + log(E(50:200)))
y = -log(material.mev) + log_K_compton;

plot(x,y)
scatter(x(50),y(50))
