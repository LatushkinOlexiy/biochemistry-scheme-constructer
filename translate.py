import argostranslate.package
import argostranslate.translate

def eng_to_ukr(text="Hello World!"):
    from_code = "en"
    to_code = "uk"

    # Download and install Argos Translate package
    argostranslate.package.update_package_index()
    available_packages = argostranslate.package.get_available_packages()
    available_package = list(
        filter(
            lambda x: x.from_code == from_code and x.to_code == to_code, available_packages
        )
    )[0]
    download_path = available_package.download()
    argostranslate.package.install_from_path(download_path)

    # Translate
    installed_languages = argostranslate.translate.get_installed_languages()
    from_lang = list(filter(
            lambda x: x.code == from_code,
            installed_languages))[0]
    to_lang = list(filter(
            lambda x: x.code == to_code,
            installed_languages))[0]
    translation = from_lang.get_translation(to_lang)
    translatedText = translation.translate(text)
    return(translatedText)
def ukr_to_eng(text="Hello World!"):
    if (text.isascii()):
        return(text)
    from_code = "uk"
    to_code = "en"

    # Download and install Argos Translate package
    argostranslate.package.update_package_index()
    available_packages = argostranslate.package.get_available_packages()
    available_package = list(
        filter(
            lambda x: x.from_code == from_code and x.to_code == to_code, available_packages
        )
    )[0]
    download_path = available_package.download()
    argostranslate.package.install_from_path(download_path)

    # Translate
    installed_languages = argostranslate.translate.get_installed_languages()
    from_lang = list(filter(
            lambda x: x.code == from_code,
            installed_languages))[0]
    to_lang = list(filter(
            lambda x: x.code == to_code,
            installed_languages))[0]
    translation = from_lang.get_translation(to_lang)
    translatedText = translation.translate(text)
    return(translatedText)

def amino_acids(text):
    text = text.lower()
    amino_acid_list = ["alanine","arginine","asparagine","aspartic acid","cysteine","glutamine","glutamate","glycine","histidine","isoleucine","leucine","lysine","methionine","phenylalanine","proline","serine","threonine","tryptophan","tyrosine","valine"]
    if text in amino_acid_list:
        return("l-"+text)
    return(text)

if __name__ == "__main__":
    print(amino_acids("alanine"))
    print(amino_acids("Methionine"))
    print(amino_acids("glucose"))
    print(ukr_to_eng("glycerine"))

