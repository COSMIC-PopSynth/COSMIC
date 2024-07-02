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
});