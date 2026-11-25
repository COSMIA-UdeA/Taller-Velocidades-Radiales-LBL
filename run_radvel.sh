SETUP=TOI233_MX_BR_ESP_1GP_ecross_RADVEL_EXAMPLE.py

radvel fit -s $SETUP
radvel plot -t rv -s $SETUP
radvel mcmc -s $SETUP
radvel derive -s $SETUP
radvel ic -t e nplanets trend -s $SETUP
radvel plot -t rv corner trend derived -s $SETUP
radvel report -s $SETUP
