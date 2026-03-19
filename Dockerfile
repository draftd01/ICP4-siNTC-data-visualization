# Use the official R base image from the Rocker project
FROM rocker/r-ver:4.3.1

# Install system dependencies for Shiny and Shiny Server
RUN apt-get update && apt-get install -y \
    sudo \
    wget \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libpng-dev \
    libtiff5-dev \
    && apt-get clean

# Install R packages needed for Shiny and Shiny Server
RUN R -e "install.packages(c('shiny', 'rmarkdown','ggplot2', 'plotly', 'heatmaply', 'bslib', 'DT', 'shinyWidgets', 'shinyjs', 'httr', 'jsonlite', 'xdvir'), repos='https://cloud.r-project.org/')"

# Download and install Shiny Server
RUN wget https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-1.5.16.958-amd64.deb && \
    gdebi -n shiny-server-1.5.16.958-amd64.deb && \
    rm shiny-server-1.5.16.958-amd64.deb

# Create a directory for the Shiny app
RUN mkdir -p /srv/shiny-server/

# Copy your custom Shiny Server config
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Expose the port for Shiny Server
EXPOSE 3839

# Run Shiny Server
CMD ["/usr/bin/shiny-server"]

