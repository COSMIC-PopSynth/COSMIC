document.addEventListener("DOMContentLoaded", function() {
    document.getElementById("hide-this-maths").classList.add("hide-me");

    document.querySelectorAll(".options-expander").forEach(expander => {
        expander.addEventListener("click", function() {
            expander.parentElement.querySelector(".options").classList.toggle("hide");
            expander.querySelector(".fa").classList.toggle("fa-chevron-down");
            expander.querySelector(".fa").classList.toggle("fa-chevron-up");
        });
    });

    document.querySelectorAll(".col-3 label").forEach(label => {
        label.addEventListener("click", function() {
            label.parentElement.querySelector("input").click();
        });
    });

    document.querySelectorAll("h3:not(.setting-group-title)").forEach(h3 => {
        let my_headings = document.querySelectorAll("h2.setting-group-title");
        for (let heading of my_headings) {
            if (heading.innerText == h3.innerText) {
                h3.style = "height: 0; margin: 0; overflow: hidden;";
                break;
            }
        }
    });

    const lis = document.querySelectorAll(".toc-drawer li");
    const bse_subheadings = document.querySelectorAll(".settings-section h2");
    for (let li of lis) {
        if (li.innerText.toLowerCase() == "binary physics") {
            new_ul = document.createElement("ul");
            for (let subheading of bse_subheadings) {
                new_li = document.createElement("li");
                new_li.innerHTML = `<a class="reference internal" href="#${subheading.id}">${subheading.innerText}</a>`;
                new_ul.appendChild(new_li);
            }
            li.appendChild(new_ul);
        }
    }

    document.querySelectorAll("input[type='checkbox']").forEach(checkbox => {
        checkbox.addEventListener("change", function() {
            const checked = checkbox.parentElement.parentElement.querySelectorAll("input[type='checkbox']");
            value = "[]";
            if (checked.length > 0) {
                value = "[";
                for (let check of checked) {
                    if (check.checked) {
                        value += `${check.value}, `;
                    }
                }
                value = value == "[" ? "[]" : value.slice(0, -2) + "]";
            }
            checkbox.parentElement.parentElement.querySelector("input.hide-me").value = value;
        });
    });

    document.querySelectorAll(".form-control").forEach(input => {
        input.addEventListener("change", function() {
            construct_files()
        });
    });

    document.getElementById("ini-file-comments").addEventListener("click", function() {
        this.classList.toggle("active");
        this.innerText = this.classList.contains("active") ? "Hide Comments" : "Show Comments";
        construct_files();
    });

    construct_files();

});


function construct_files() {
    let ini_file = "<span class='c1'>; COSMIC INI file\n</span>"
    let BSE_dict = ""
    let BSE_dict_string = "<span class='p'>BSE_settings </span><span class='o'> = </span><span class='p'>{</span>"

    reached_bse = false;

    const include_comments = document.getElementById("ini-file-comments").classList.contains("active");

    const els = document.querySelectorAll(".setting-card, .settings-section, .form-control");
    for (let el of els) {
        if (el.classList.contains("setting-card")) {
            if (el.getAttribute("data-category") == "bse") {
                reached_bse = true;
            }
            ini_file += `<span class='k'>\n[${el.getAttribute("data-category")}]\n</span>`;
        } else if (el.classList.contains("settings-section")) {
            subheading = el.querySelector("h2").innerText;
            ini_file += '<span class="c1">\n;;;;' + ';'.repeat(subheading.length) + ';;;;\n</span>';
            ini_file += `<span class="c1">;;; ${subheading} ;;;\n</span>`;
            ini_file += '<span class="c1">;;;;' + ';'.repeat(subheading.length) + ';;;;\n</span>';
        } else if (el.classList.contains("form-control")) {
            const full_setting = el.parentElement.parentElement.parentElement;
            if (include_comments) {
                ini_file += `<span class="c1">\n; ${full_setting.querySelector(".name").innerText}\n</span>`;

                // include the description as a comment
                description = full_setting.querySelector(".description").innerText.replaceAll('\n', '');
                ini_file += `<span class="c1">; ${description}\n</span>`;

                // include the options-preface as a comment
                options_preface = full_setting.querySelector(".options-preface");
                if (options_preface && options_preface.innerText != "") {
                    ini_file += `<span class='c1'>; ${options_preface.innerText.replaceAll('\n', '')}\n</span>`;
                }

                // include the options as a comment
                options = full_setting.querySelectorAll("li");
                if (options.length > 0) {
                    ini_file += "<span class='c1'>; Options: \n";
                    for (let option of options) {
                        ini_file += `<span class='c1'>;    ${option.innerText.replaceAll('\n', '')}\n`;
                    }
                    ini_file += "</span>";
                }

                // include the default value as a comment
                default_value = full_setting.querySelector(".default").innerText;
                ini_file += `<span class="c1">; ${default_value}\n</span>`;
            }
            label = full_setting.querySelector(".name").innerText;
            value = el.value;
            ini_file += `<span class="na">${label}</span><span class="na"> = </span><span class="s">${value}\n</span>`;
            if (reached_bse) {
                BSE_dict[label] = value;
                BSE_dict_string += `<span class="s1">'${label}'</span><span class="p">: </span><span class="mf">${value}</span><span class="p">, </span>`
            }
        }
    }

    document.querySelector("#ini-file .highlight-ini pre").innerHTML = ini_file;
    document.querySelector("#python-bse-settings-dictionary .highlight-python pre").innerHTML = BSE_dict_string.slice(0, -9) + "</span><span class='p'>}</span>";
}
