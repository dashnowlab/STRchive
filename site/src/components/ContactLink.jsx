import Button from "@/components/Button";
import ContactForm from "@/components/ContactForm";
import Dialog from "@/components/Dialog";
import classes from "./Contact.module.css";

const ContactLink = () => (
  <Dialog
    title="Contact Us"
    trigger={
      <Button
        id="footer-contact"
        className={classes.link}
        data-tooltip="Contact us"
      >
        Contact
      </Button>
    }
  >
    <ContactForm />
  </Dialog>
);

export default ContactLink;
