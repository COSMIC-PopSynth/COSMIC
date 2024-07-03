document.addEventListener("DOMContentLoaded", function() {
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

    construct_files()
});


function construct_files() {
    let ini_file = "<span class='c1'>; COSMIC INI file\n</span>"
    let BSE_dict = {}

    const els = document.querySelectorAll(".setting-card, .settings-section, .form-control");
    for (let el of els) {
        if (el.classList.contains("setting-card")) {
            ini_file += `<span class='k'>\n[${el.getAttribute("data-category")}]\n</span>`;
        } else if (el.classList.contains("settings-section")) {
            subheading = el.querySelector("h2").innerText;
            ini_file += '<span class="c1">\n;;;;' + ';'.repeat(subheading.length) + ';;;;\n</span>';
            ini_file += `<span class="c1">;;; ${subheading} ;;;\n</span>`;
            ini_file += '<span class="c1">;;;;' + ';'.repeat(subheading.length) + ';;;;\n</span>';
        } else if (el.classList.contains("form-control")) {
            label = el.parentElement.parentElement.querySelector(".name").innerText;
            value = el.value;
            ini_file += `<span class="na">${label}</span><span class="na"> = </span><span class="s">${value}\n</span>`;
        }
    }

    document.querySelector("#ini-file .highlight-ini pre").innerHTML = ini_file;
}