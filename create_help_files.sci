//Criar os arquivos .xml das funções em path+'help_xmls'
path=get_absolute_file_path()
//disp(path)
chdir(path+'macros')
mkdir(path+'help')
help_from_sci(path+'macros',path+'help')

//Cria visualização do help em path+'help'
chdir(path)
mkdir('help')
chdir(path+'help')
xmltojar(path+'help','Internal Fluid Flow','en_US')

//Cria arquivo path+'jar\scilab_en_US_help.jar'
//chdir(path)
//mkdir('jar')
//jcreatejar(path+'jar\scilab_en_US_help.jar',..
//           path+'help')

add_help_chapter("Internal Fluid Flow",path)//+'jar')
