FROM bgruening/galaxy-stable

ENV GALAXY_CONFIG_BRAND LungCARD

WORKDIR /galaxy-central

# Add galaxy.yml with ports configured and admin_user
ADD galaxy.yml /galaxy-central/config/galaxy.yml
ADD dependency_resolvers_conf.xml /galaxy-central/config/dependency_resolvers_conf.xml

RUN add-tool-shed --url 'http://testtoolshed.g2.bx.psu.edu/' --name 'Test Tool Shed'

# Adding the tool definitions to the container
ADD my_tool_list.yml $GALAXY_ROOT/my_tool_list.yml

# Install my_tools_list
RUN install-tools $GALAXY_ROOT/my_tool_list.yml

# Add lungCARD tools

ADD /tools/add_new_variant /galaxy-central/tools/add_new_variant 
ADD /tools/report /galaxy-central/tools/report
ADD /tools/snpEff /galaxy-central/tools/snpEff
ADD /tools/vardict /galaxy-central/tools/vardict

# Add bed_files and .loc files (with container paths)

ADD /bed_files /galaxy-central/tool-data/bed_files
ADD bed_files.loc /galaxy-central/tool-data/bed_files.loc
ADD snpEff_genomes.loc /galaxy-central/tool-data/snpEff_genomes.loc

# Add tool_data_table_conf.xml and tool_conf.xml

ADD /config/tool_data_table_conf.xml /galaxy-central/config/tool_data_table_conf.xml

ADD /config/tool_conf.xml /galaxy-central/config/tool_conf.xml

# Add workfolws

ADD /workflows /galaxy-central/workflows

# Add static files

ADD /static/welcome.html /galaxy-central/static/welcome.html
ADD /static/logo.png /galaxy-central/static/logo.png
ADD /static/welcome.html /etc/galaxy/web/welcome.html
ADD /static/logo.png /etc/galaxy/web/logo.png


